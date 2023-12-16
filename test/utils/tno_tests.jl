@testset "TNO Yellow Book Equation Sets" begin

    @testset "Pasquill-Gifford dispersion tests" begin

        # Plume dispersion
        known = [(ClassA, 0.6170244776546403, 0.3299295029728027),
                 (ClassB, 0.43445506897745934, 0.2685541668129271),
                 (ClassC, 0.24613414407041376, 0.2545468210566936),
                 (ClassD, 0.15096247150079659, 0.2297247262425074),
                 (ClassE, 0.11551744100090941, 0.1713537324266559),
                 (ClassF, 0.07661871086795012, 0.135591567342651)]
        @testset "Plume dispersion, stability class $class" for (class,sy,sz) in known
            @test crosswind_dispersion(1.2, Plume, class, TNO()) ≈ sy
            @test vertical_dispersion(1.2, Plume, class, TNO()) ≈ sz
        end

        # Puff dispersion
        known = [(ClassA, 0.30851223882732015, 0.336, 0.156),
                 (ClassB, 0.21722753448872967, 0.276, 0.156),
                 (ClassC, 0.12306707203520688, 0.264, 0.156),
                 (ClassD, 0.07548123575039829, 0.240, 0.156),
                 (ClassE, 0.057758720500454705, 0.18, 0.156),
                 (ClassF, 0.03830935543397506, 0.144, 0.156)]
       @testset "Puff dispersion, stability class $class" for (class,sy,sz,sx) in known
            @test downwind_dispersion(1.2, Puff, class, TNO()) ≈ sx
            @test crosswind_dispersion(1.2, Puff, class, TNO()) ≈ sy
            @test vertical_dispersion(1.2, Puff, class, TNO()) ≈ sz
       end
    end

    @testset "Windspeed (default)" begin

        u0, z0 = 3.0, 1.0
        a = SimpleAtmosphere(windspeed=u0, windspeed_height=z0, stability=ClassA)
        s = Scenario(Substance(:null,0,0,0,0,0,0,0,0),HorizontalJet(0,0,0,0,1.0,0,0,0),a)
        @test _windspeed(s,10,TNO()) == _windspeed(a,10,TNO()) == _windspeed(u0,z0,10,ClassA,TNO())
        @test _windspeed(u0,z0,10,ClassA,TNO()) == _windspeed(u0,z0,10,ClassA,DefaultSet())
    end
end