@testset "SLAB puff tests" begin

    @test GasDispersion._slab_stab.([ClassA,ClassB,ClassC,ClassD,ClassE,ClassF]) ≈ [1.0,2.0,3.0,4.0,5.0,6.0]

@testset "Antoine Coefficent Recovery" begin
    Base.isapprox(a::GasDispersion.Antoine, b::GasDispersion.Antoine) = all([
        getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
        if typeof(getproperty(a,k))<:Number ])

    s = Substance(name=:test,molar_weight=1,vapor_pressure=1,gas_density=1,
                  liquid_density=1,boiling_temp =1,latent_heat=1,gas_heat_capacity=1,
                  liquid_heat_capacity=1) 
    @test GasDispersion._slab_antoine(s) ≈ GasDispersion.Antoine(0.0,-1.0,0.0)

    # propane from Perry's 8th edition, DIPPR correlations
    pv = GasDispersion.DIPPRVaporPressure(59.078,-3_492.6,-6.0669,1.0919e-5,2)
    ρl = GasDispersion.DIPPRLiquidDensity(0.044096,369.83,1.3757,0.27453,369.83,0.29359)
    Δhv = GasDispersion.DIPPRLatentHeat(0.044096,369.83,2.9209e7,0.78237,-0.77319,0.39246,0)
    cp_ig = GasDispersion.DIPPRIdealGasHeatCapacity(0.044096,369.83,0.5192e5,1.9245e5,1.6265e3,1.168e5,723.6)
    cp_l = GasDispersion.DIPPRLiquidHeatCapacity(GasDispersion.Eq2,0.044096,369.83,62.983,113_630,633.21,-873.46,0)

    propane = Substance(name="propane",
                        molar_weight=0.044096,
                        vapor_pressure=pv,
                        gas_density=nothing,
                        liquid_density=ρl,
                        reference_temp=288.15,
                        reference_pressure=101325,
                        k=1.3,
                        boiling_temp=231.02,
                        latent_heat=Δhv,
                        gas_heat_capacity=cp_ig,
                        liquid_heat_capacity=cp_l)
    @test GasDispersion._slab_antoine(propane) ≈ GasDispersion.Antoine(20.75136471777637, 1929.0514351993063, -21.974809894382986)

    propane = Substance(name="propane",
                        molar_weight=0.044096,
                        vapor_pressure=nothing,
                        gas_density=nothing,
                        liquid_density=ρl,
                        reference_temp=288.15,
                        reference_pressure=101325,
                        k=1.3,
                        boiling_temp=231.02,
                        latent_heat=Δhv,
                        gas_heat_capacity=cp_ig,
                        liquid_heat_capacity=cp_l)
    @test GasDispersion._slab_antoine(propane) ≈ GasDispersion.Antoine(8.086226333190652, 1868.0800074937044, 0.0)

end

@testset "INPR2 Horizontal Jet" begin
# this directly tests the slab.jl submodule against the given test problem 2
    inp = GasDispersion.SLAB_Input(idspl =  2,
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
    res = GasDispersion.slab.slab_main(inp)
    out = readdlm("test_data/slab_inpr2_out.txt", Float64; header=false);
  
    # @testset "release gas properties" begin
    #     @test res.p.rgp.wms   ≈ 1.7031E-02
    #     @test res.p.rgp.cps   ≈ 2.0459E+03
    #     @test res.p.rgp.ts    ≈ 2.3957E+02
    #     @test res.p.rgp.rhos  ≈ 0.8663594467626095
    #     @test res.p.rgp.tbp   ≈ 2.3957E+02
    #     @test res.p.rgp.cmed0 ≈ 8.1000E-01
    #     @test res.p.rgp.cpsl  ≈ 4.6118E+03
    #     @test res.p.rgp.dhe   ≈ 1.1700E+06
    #     @test res.p.rgp.rhosl ≈ 6.0300E+02
    #     @test res.p.rgp.spa   ≈ 12.422298284426265
    #     @test res.p.rgp.spb   ≈ 2976.01
    #     @test res.p.rgp.spc   ≈ 0.0000E+00
    # end


    # @testset "spill characteristics" begin
    #     @test res.p.spl.idspl≈ 2
    #     @test res.p.spl.qs   ≈ 1.0787E+02
    #     @test res.p.spl.tsd  ≈ 3.8100E+02
    #     @test res.p.spl.qtcs ≈ 41098.47
    #     @test res.p.spl.qtis ≈ 0.0000E+00
    #     @test res.p.spl.as   ≈ 9.3000E-01
    #     @test res.p.spl.ws   ≈ 0.0000E+00
    #     @test res.p.spl.bs   ≈ 0.4821825380496478
    #     @test res.p.spl.hs   ≈ 1.0000E+00
    #     @test res.p.spl.us   ≈ 25.59323553673247
    # end


    # @testset "field parameters" begin
    #     @test res.p.fld.tav  ≈  1.0000E+01
    #     @test res.p.fld.hmx  ≈  726.038044311263
    #     @test res.p.fld.xffm ≈  2.8000E+03
    #     @test res.p.fld.zp ==  [0.0, 1.0, 0.0, 0.0]
    # end


    # @testset "ambient meteorological properties" begin
    #     @test res.p.met.wmae  ≈  0.0288353038537132
    #     @test res.p.met.cpaa  ≈  1011.8542876036452
    #     @test res.p.met.rhoa  ≈  1.147650750527487
    #     @test res.p.met.za    ≈  2.0000E+00
    #     @test res.p.met.pa    ≈  101325.0
    #     @test res.p.met.ua    ≈  4.5000E+00
    #     @test res.p.met.ta    ≈  3.0620E+02
    #     @test res.p.met.rh    ≈  2.1300E+01
    #     @test res.p.met.uastr ≈  0.26631640297328646
    #     @test res.p.met.stab  ≈  4.518466475995308
    #     @test res.p.met.ala   ≈  2.2100E-02
    #     @test res.p.met.z0    ≈  3.0000E-03
    # end

    # @testset "additional parameters" begin
    #     @test res.p.xtra.ncalc == 1
    #     @test res.p.xtra.nssm  == 3
    #     @test res.p.xtra.grav ≈ 9.80665
    #     @test res.p.xtra.rr ≈ 8.31431
    #     @test res.p.xtra.xk ≈ 0.41

    # end

    # instantaneous spatially averaged cloud parameters
    result = [res.s.x res.s.zc res.s.h res.s.bb res.s.b res.s.bbx res.s.bx res.s.cv res.s.rho res.s.t res.s.u res.s.uab res.s.cm res.s.cmev res.s.cmda res.s.cmw res.s.cmwv res.s.wc res.s.vg res.s.ug res.s.w res.s.v res.s.vx ]
    @test result ≈ out[:, 1:23]

    #averaged cloud parameters
    result = [res.cc.cc res.cc.betac res.cc.sig res.cc.t res.cc.xc res.cc.betax res.cc.tim res.cc.tcld res.cc.bbc]
    @test result ≈ out[:,24:32]

end

@testset "Burro LNG test" begin
# this tests the output against an example from the Burro LNG dispersion tests
# the output of slab.jl was compared to that generated from SLAB (fortran) using
# the same input file
    R = 8.31446261815324
    ρ = 1.76
    T = 273.15-162
    P = 101325
    MW = ρ*R*T/P

    s = Substance(name = :BurroLNG,
                molar_weight = MW,
                vapor_pressure = nothing,
                gas_density = ρ,
                liquid_density = 425.6,
                reference_temp=T,
                reference_pressure=P,
                k=1.4,
                boiling_temp = 111.66, # K, Methane,
                latent_heat = 509880,  # J/kg, Methane
                gas_heat_capacity = 2240, # J/kg/K, Methane
                liquid_heat_capacity = 3349) # J/kg/K, Methane
    r = HorizontalJet( mass_rate = (0.23*425.6),
                duration = 174,
                diameter = 1,
                velocity = 70.815,
                height = 0,
                pressure = P,
                temperature = T,
                fraction_liquid = 0)
    a = SimpleAtmosphere(windspeed=10.9, temperature=298, stability=ClassF)
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