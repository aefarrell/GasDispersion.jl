import GasDispersion.slab: SLAB_Release_Parameters, SLAB_Plume_Integrated, SLAB_Thermo_State, SLAB_Plume_State
import GasDispersion.slab: SLAB_Release_Source, SLAB_Release_Ambient_Air, SLAB_Antoine
import GasDispersion.slab: _slab_psat
import GasDispersion.slab: plume_thermodynamics
import GasDispersion.slab: total_mass_fraction, vapour_mass_fraction, water_mass_fraction, water_vapour_mass_fraction, dry_air_mass_fraction, volume_fraction
import GasDispersion.slab: total_specific_heat, transported_energy
import GasDispersion.slab: vapour_partial_pressure, vapour_check_water, vapour_check_substance, vapour_phase_initial_guess

include("../testutils/trace.jl")

const spaw = 15.08
const spbw = 5514.
const wma = 0.02896
const cpa = 1005.87
const grav = 9.80665
const pa = 101325.

const wmw = 0.01802
const rhowl = 1000.
const cpwv = 1846.
const cpwl = 4178.

@testset "SLAB thermodynamics model tests" begin
    # using the INPR2 test data to examine the thermo model
    results = readtrace("test_data/checkpoint_thermo.jls")

    calls = []
    checkpoints = Dict()
    for result in results
        if result.label == :call
            #new function call, cache states
            push!(calls,checkpoints)
            checkpoints = Dict()
        end
        checkpoints[result.label] = result.state
    end
    push!(calls,checkpoints)
    popfirst!(calls)

    for call in filter((c)-> c[:call].idpf !=2, calls)
        rls_params,idpf,xn,timn,rmi,t0,cmev0,cm0,cmw0,cmwv0,cp0,r,r0,sft,bse = call[:call]
        bbv,bv,g,gw,zc,qint,sfu,sfy,sfz = zeros(eltype(r),9)
        bbv0,bv0,g0,gw0,zc0,qint0,sft0,sfu0,sfy0,sfz0 = zeros(eltype(r),10)

        qs = rls_params.spl.qs
        wmae = rls_params.met.wmae
        wms = rls_params.rgp.wms
        cmwa = rls_params.met.cmwa
        cpaa = rls_params.met.cpaa
        ta = rls_params.met.ta
        cps = rls_params.rgp.cps
        cpsl = rls_params.rgp.cpsl
        ts = rls_params.rgp.ts
        tbp = rls_params.rgp.tbp
        dhe = rls_params.rgp.dhe

        release_antoine = SLAB_Antoine(rls_params.rgp.spa, rls_params.rgp.spb, rls_params.rgp.spc)
        water_antoine = SLAB_Antoine(spaw, spbw, 0.0)
        source = SLAB_Release_Source(bse, qs, wms, cps, cpsl, ts, tbp, dhe, release_antoine)
        ambient_air = SLAB_Release_Ambient_Air(wmae, wma, wmw, cmwa, cpaa, cpa, cpwv, cpwl, ta, water_antoine)
        p = SLAB_Release_Parameters(source, ambient_air)
        curr = SLAB_Plume_Integrated(r,bbv,bv,g,gw,zc,qint,sft,sfu,sfy,sfz)
        last = SLAB_Plume_State(
            SLAB_Plume_Integrated(r0,bbv0,bv0,g0,gw0,zc0,qint0,sft0,sfu0,sfy0,sfz0),
            SLAB_Thermo_State(cm0,0.0,cmw0,1.0-cm0-cmw0,cmwv0,cmw0-cmwv0,cmev0,cm0-cmev0,cp0,t0,0.0,nothing,0.0)
        )

        new_thermo = plume_thermodynamics(xn, curr, last, p)
    
        # testing concentrations
        concentrations = call[:concentrations]
        cm = new_thermo.cm
        @test cm ≈ concentrations.cm
        cv = new_thermo.cv
        @test cv ≈ concentrations.cv
        cmw = new_thermo.cmw
        @test cmw ≈ concentrations.cmw
        cmda = new_thermo.cmda
        @test cmda ≈ concentrations.cmda

        cmwv = water_vapour_mass_fraction(new_thermo.cmw, curr, last) 
        @test cmwv ≈ concentrations.cmwv
        cmwd = cmw - cmwv
        @test cmwd ≈ concentrations.cmwd
        
        cmev = vapour_mass_fraction(new_thermo.cm, curr, last)
        @test cmev ≈ concentrations.cmev
        cmed = cm - cmev
        @test cmed ≈ concentrations.cmed

        # testing transported energy and initial thermodynamic state
        transported_energy = call[:transported_energy]
        etrn = new_thermo.transported_energy
        @test etrn ≈ transported_energy.etrn
        cp = total_specific_heat(cmda, cmwv, cmwd, cmev, cmed, p)
        @test cp ≈ transported_energy._cp
        T = etrn/cp
        @test T ≈ transported_energy.t
        
        # testing initial phase check
        initial_phase_check = call[:initial_phase_check]
        sub = p.source.antoine
        @test _slab_psat(T, sub) ≈ initial_phase_check.fft
        wtr = p.ambient_air.antoine
        @test _slab_psat(T, wtr) ≈ initial_phase_check.ggt

        @test vapour_partial_pressure(cm, wms, cmev, cmwv, cmda, p) ≈ initial_phase_check.pe
        @test vapour_partial_pressure(cmw, wmw, cmev, cmwv, cmda, p) ≈ initial_phase_check.pw

        @test vapour_check_substance(cm, cmev, cmwv, cmda, T, last, p) == (initial_phase_check.ie == 1)
        @test vapour_check_water(cmw, cmw, cmev, cmda, T, last, p) == (initial_phase_check.iw == 1)

        # testing vapour phase change initialization
        if :vapour_temperature_initial in keys(call)
            vapour_temperature_initial = call[:vapour_temperature_initial]
            dlti = vapour_temperature_initial.dlti_
            dltm = vapour_temperature_initial.dltm
            @test max(dltm,dlti) ≈ vapour_phase_initial_guess(cm,cmev,cmwv,cmda,T,cp,p)
        end
    end

end