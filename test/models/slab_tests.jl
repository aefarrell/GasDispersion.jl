include("../../src/models/slab/slab.jl")

using .slab

# this is incredibly verbose to start, I want to make sure it is reproducing the 
# SLAB output file entirely later I can probably pare this back to just ensure 
# every branch is being hit

@testset "SLAB Horizontal Jet INPR2" begin
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

    # the following checks the horizontal jet source
    # initialization (label 400)
    @testset "horizontal jet initialization" begin
        @test res.ic.x[1]    ≈ 1
        @test res.ic.vg[1]   ≈ 0
        @test res.ic.wc[1]   ≈ 0
        @test res.ic.ug[1]   ≈ 0
        @test res.ic.bb[1]   ≈ 0.482182533 
        @test res.ic.b[1]    ≈ 0.433964282 
        @test res.ic.rho[1]  ≈ 4.532027501773091
        @test res.ic.h[1]    ≈ 0.964365065
        @test res.ic.uab[1]   ≈ 3.4753698069472003
        @test res.ic.u[1]    ≈ 25.59323553673247
        @test res.ic.w[1]    ≈ 0.6154294519602936
        @test res.ic.v[1]    ≈ 0.49637768050087083
        @test res.ic.vx[1]   ≈ 0
    end

    # the following checks the near field region 
    # calculation (label 600)
    @testset "steady state plume dispersion, near field calculation" begin
        @test res.ic.x[54]  ≈ 848.5204149435892
        @test res.ic.cm[54] ≈ 0.0037592610280749433
    end

end