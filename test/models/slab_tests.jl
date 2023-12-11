include("../../src/models/slab/slab.jl")

using .slab
using DelimitedFiles: readdlm

@testset "SLAB puff tests" begin

    @test GasDispersion._slab_stab.([ClassA,ClassB,ClassC,ClassD,ClassE,ClassF]) ≈ [1.0,2.0,3.0,4.0,5.0,6.0]

@testset "INPR2 Horizontal Jet" begin
# this directly tests the slab.jl submodule against the given test problem 2
    inp = SLAB_Input(idspl =  2,
                     ncalc =  1,
                     wms   =  0.017031,
                     cps   =  2045.90,
                     tbp   =  239.57,
                     cmed0 =  0.81,
                     dhe   =  1170000.0,
                     cpsl  =  4611.80,
                     rhosl =  603.00,
                     spb   =  2976.01,
                     spc   =  0.00,
                     ts    =  239.57,
                     qs    =  107.87,
                     as    =  0.93,
                     tsd   =  381.0,
                     qtis  =  0.00,
                     hs    =  1.00,
                     tav   =  10.00,
                     xffm  =  2800.00,
                     zp    =  [0.00, 1.00, 0.00, 0.00],
                     z0    =  0.003000,
                     za    =  2.00,
                     ua    =  4.50,
                     ta    =  306.20,
                     rh    =  21.30,
                     stab  =  0.00,
                     ala   =  0.0221)
    res = slab_main(inp)
    out = readdlm("test_data/slab_inpr2_out.txt", Float64; header=false);
    # These test whether or not the parameters calculated
    # from SLAB.for up to line 567 (call editin)
    # are doing what they are supposed to be doing
    @testset "release gas properties" begin
        @test res.p.rgp.wms   ≈ 1.7031E-02
        @test res.p.rgp.cps   ≈ 2.0459E+03
        @test res.p.rgp.ts    ≈ 2.3957E+02
        @test res.p.rgp.rhos  ≈ 0.8663594467626095
        @test res.p.rgp.tbp   ≈ 2.3957E+02
        @test res.p.rgp.cmed0 ≈ 8.1000E-01
        @test res.p.rgp.cpsl  ≈ 4.6118E+03
        @test res.p.rgp.dhe   ≈ 1.1700E+06
        @test res.p.rgp.rhosl ≈ 6.0300E+02
        @test res.p.rgp.spa   ≈ 12.422298284426265
        @test res.p.rgp.spb   ≈ 2976.01
        @test res.p.rgp.spc   ≈ 0.0000E+00
    end


    @testset "spill characteristics" begin
        @test res.p.spl.idspl≈ 2
        @test res.p.spl.qs   ≈ 1.0787E+02
        @test res.p.spl.tsd  ≈ 3.8100E+02
        @test res.p.spl.qtcs ≈ 41098.47
        @test res.p.spl.qtis ≈ 0.0000E+00
        @test res.p.spl.as   ≈ 9.3000E-01
        @test res.p.spl.ws   ≈ 0.0000E+00
        @test res.p.spl.bs   ≈ 0.4821825380496478
        @test res.p.spl.hs   ≈ 1.0000E+00
        @test res.p.spl.us   ≈ 25.59323553673247
    end


    @testset "field parameters" begin
        @test res.p.fld.tav  ≈  1.0000E+01
        @test res.p.fld.hmx  ≈  726.038044311263
        @test res.p.fld.xffm ≈  2.8000E+03
        @test res.p.fld.zp ==  [0.0, 1.0, 0.0, 0.0]
    end


    @testset "ambient meteorological properties" begin
        @test res.p.met.wmae  ≈  0.0288353038537132
        @test res.p.met.cpaa  ≈  1011.8542876036452
        @test res.p.met.rhoa  ≈  1.147650750527487
        @test res.p.met.za    ≈  2.0000E+00
        @test res.p.met.pa    ≈  101325.0
        @test res.p.met.ua    ≈  4.5000E+00
        @test res.p.met.ta    ≈  3.0620E+02
        @test res.p.met.rh    ≈  2.1300E+01
        @test res.p.met.uastr ≈  0.26631640297328646
        @test res.p.met.stab  ≈  4.518466475995308
        @test res.p.met.ala   ≈  2.2100E-02
        @test res.p.met.z0    ≈  3.0000E-03
    end

    @testset "additional parameters" begin
        @test res.p.xtra.ncalc == 1
        @test res.p.xtra.nssm  == 3
        @test res.p.xtra.grav ≈ 9.80665
        @test res.p.xtra.rr ≈ 8.31431
        @test res.p.xtra.xk ≈ 0.41

    end

    @testset "instantaneous spatially averaged cloud parameters" begin
        @test res.s.x ≈ out[:, 1];
        @test res.s.zc ≈ out[:, 2];
        @test res.s.h ≈ out[:, 3];
        @test res.s.bb ≈ out[:, 4];
        @test res.s.b ≈ out[:, 5];
        @test res.s.bbx ≈ out[:, 6];
        @test res.s.bx ≈ out[:, 7];
        @test res.s.cv ≈ out[:, 8];
        @test res.s.rho ≈ out[:, 9];
        @test res.s.t ≈ out[:, 10];
        @test res.s.u ≈ out[:, 11];
        @test res.s.uab ≈ out[:, 12];
        @test res.s.cm ≈ out[:, 13];
        @test res.s.cmev ≈ out[:, 14];
        @test res.s.cmda ≈ out[:, 15];
        @test res.s.cmw ≈ out[:, 16];
        @test res.s.cmwv ≈ out[:, 17];
        @test res.s.wc ≈ out[:, 18];
        @test res.s.vg ≈ out[:, 19];
        @test res.s.ug ≈ out[:, 20];
        @test res.s.w ≈ out[:, 21];
        @test res.s.v ≈ out[:, 22];
        @test res.s.vx ≈ out[:, 23];
    end

    @testset "time averaged cloud parameters" begin
        @test res.cc.x ≈ out[:, 1];
        @test res.cc.cc ≈ out[:, 24];
        @test res.cc.b ≈ out[:, 5];
        @test res.cc.betac ≈ out[:, 25];
        @test res.cc.zc ≈ out[:, 2];
        @test res.cc.sig ≈ out[:, 26];
        @test res.cc.t ≈ out[:, 27];
        @test res.cc.xc ≈ out[:, 28];
        @test res.cc.bx ≈ out[:, 7];
        @test res.cc.betax ≈ out[:, 29];
        @test res.cc.tim ≈ out[:, 30];
        @test res.cc.tcld ≈ out[:, 31];
        @test res.cc.bbc ≈ out[:, 32];
    end
end

@testset "Burro LNG test" begin
# this tests the output against an example from the Burro LNG dispersion tests
# the output of slab.jl was compared to that generated from SLAB (fortran) using
# the same input file

    s = Substance(name = :BurroLNG,
                gas_density = 1.76,
                liquid_density = 425.6,
                reference_temp=(273.15-162),
                reference_pressure=101325.0,
                boiling_temp = 111.66, # K, Methane,
                latent_heat = 509880.0,  # J/kg, Methane
                gas_heat_capacity = 2240.0, # J/kg/K, Methane
                liquid_heat_capacity = 3349.0) # J/kg/K, Methane
    r = Release( mass_rate = (0.23*425.6),
                duration = 174,
                diameter = 1.0,
                velocity = 70.815,
                height = 0.0,
                pressure = 101325.0,
                temperature = (273.15-162),
                fraction_liquid = 0.0)
    a = DryAir(windspeed=10.9, temperature=298, stability=ClassF)
    scn = Scenario(s,r,a)
    rls = puff(scn, SLAB; )

    # basic checks that it returns the expected object
    @test isa(rls,GasDispersion.SLABSolution)
    @test isa(rls, Puff)

    # basic domain checks
    @test rls(0.0,0.0,r.h,1.0) ≈ 1.0
    @test rls(-1.0,0.0,0.0,1.0) == 0.0
    @test_throws ErrorException rls(rls.in.xffm + eps(Float64), 0.0, 0.0, 1.0)

    # check against time averaged cloud concentrations
    burro = readdlm("test_data/slab_burro_z0.txt", Float64; header=false);
    x = burro[:,1]
    t = burro[:,2]
    tcld = burro[:,3]
    bbc = burro[:,4]
    cs = burro[:,5:end]
    ms = [0.0,0.5,1.0,1.5,2.0,2.5]
    ys = bbc*ms'
    result = hcat([[ rls(x[j],ys[j,i],0.0,t[j]) for j in 1:60 ] for i in 1:6]...)
    @test result ≈ cs
    
end

end