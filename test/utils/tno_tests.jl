@testset "TNO Yellow Book Equation Sets" begin

    @testset "Windspeed tests" begin
        u0, z0 = 2.0, 10.0
        known = [(ClassA(), 15.0, 2.105802252736218),
                 (ClassB(), 15.0, 2.1134438908570603),
                 (ClassC(), 15.0, 2.1298117227821707),
                 (ClassD(), 15.0, 2.176091259055681),
                 (ClassE(), 15.0, 2.3401491268457204),
                 (ClassF(), 15.0, 2.542796254727042)]
        @testset "Windspeed, stability class $class" for (class, za, ua) in known
            atm = SimpleAtmosphere(;windspeed=u0,windspeed_height=z0,stability=class)
            @test GasDispersion.windspeed(atm, za, TNOPlume) ≈ ua
        end
    end

    @testset "Pasquill-Gifford dispersion tests" begin

        # Plume dispersion
        known = [(ClassA(), 1.2, 0.6170244776546403, 0.3299295029728027),
                 (ClassB(), 1.2, 0.43445506897745934, 0.2685541668129271),
                 (ClassC(), 1.2, 0.24613414407041376, 0.2545468210566936),
                 (ClassD(), 1.2, 0.15096247150079659, 0.2297247262425074),
                 (ClassE(), 1.2, 0.11551744100090941, 0.1713537324266559),
                 (ClassF(), 1.2, 0.07661871086795012, 0.135591567342651)]
        @testset "Plume dispersion, stability class $class" for (class,x,sy,sz) in known
            @test GasDispersion.crosswind_dispersion(x, class, TNOPlume) ≈ sy
            @test GasDispersion.vertical_dispersion(x, class, TNOPlume) ≈ sz
        end

        # Puff dispersion
        known = [(ClassA(), 1.2, 0.30851223882732015, 0.336, 0.156),
                 (ClassB(), 1.2, 0.21722753448872967, 0.276, 0.156),
                 (ClassC(), 1.2, 0.12306707203520688, 0.264, 0.156),
                 (ClassD(), 1.2, 0.07548123575039829, 0.240, 0.156),
                 (ClassE(), 1.2, 0.057758720500454705, 0.18, 0.156),
                 (ClassF(), 1.2, 0.03830935543397506, 0.144, 0.156)]
       @testset "Puff dispersion, stability class $class" for (class,x,sy,sz,sx) in known
            @test GasDispersion.downwind_dispersion(x, class, TNOPuff) ≈ sx
            @test GasDispersion.crosswind_dispersion(x, class, TNOPuff) ≈ sy
            @test GasDispersion.vertical_dispersion(x, class, TNOPuff) ≈ sz
       end
    end

end