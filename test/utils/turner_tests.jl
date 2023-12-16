@testset "Turner 1970 Equation Sets" begin

    @testset "Pasquill-Gifford dispersion tests" begin

        # Plume dispersion
        known = [(ClassA, 0.5787971924008366, 0.10632089486974931, 16.85072624954395, 130.3205905733347, 1.8113400926195988e10),
                 (ClassB, 0.39564838506913175, 0.1605299063149382, 12.751343710359913, 52.9324088164219, 226360.1640449754),
                 (ClassC, 0.22977233242748757, 0.13221263370966554, 8.735202083580425, 33.17271865454032, 32301.152835018),
                 (ClassD, 0.1508249156447098, 0.10858929353740095, 5.44235676295498, 18.33607966620897, 1158.7773561551285),
                 (ClassE, 0.10742276488910325, 0.09522313470055986, 4.156640589233231, 13.43619462405115, 346.73685045253166),
                 (ClassF, 0.07894741678277778, 0.0659507672737797, 2.625547335519663, 8.572359097689173, 106.16955571987263)]
        @testset "Plume dispersion, stability class $class" for (class,sy,sz1,sz2,sz3,sz4) in known
            @test crosswind_dispersion(1.2, Plume, class, Turner) ≈ sy
            @test vertical_dispersion(1.2, Plume, class, Turner) ≈ sz1
            @test vertical_dispersion(120, Plume, class, Turner) ≈ sz2
            @test vertical_dispersion(520, Plume, class, Turner) ≈ sz3
            @test vertical_dispersion(1e6, Plume, class, Turner) ≈ sz4
        end

    end

end