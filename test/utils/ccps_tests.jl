@testset "CCPS Equation Sets" begin

    @testset "Pasquill-Gifford dispersion tests" begin

        # Plume dispersion
        urban = [(ClassA, 0.38390787316433456, 0.2878273553646431),
                 (ClassB, 0.38390787316433456, 0.2878273553646431),
                 (ClassC, 0.26393666280048, 0.24),
                 (ClassD, 0.19195393658216728, 0.16769841403825084),
                 (ClassE, 0.13196833140024, 0.09591371646531512),
                 (ClassF, 0.13196833140024, 0.09591371646531512)]
        @testset "Plume dispersion, urban terrain, stability class $class" for (class,cwind,vert) in urban
            @test crosswind_dispersion(1.2, Plume, class, CCPSUrban()) ≈ cwind
            @test vertical_dispersion(1.2, Plume, class, CCPSUrban()) ≈ vert
        end

        rural = [(ClassA, 0.2639841614254575, 0.24),
                 (ClassB, 0.19198848103669636, 0.144),
                 (ClassC, 0.13199208071272875, 0.09598848207318536),
                 (ClassD, 0.09599424051834818, 0.07193528734898633),
                 (ClassE, 0.07199568038876113, 0.03598704466392099),
                 (ClassF, 0.04799712025917409, 0.019193090487424527)]
        @testset "Plume dispersion, rural terrain, stability class $class" for (class,cwind,vert) in rural
            @test crosswind_dispersion(1.2, Plume, class, CCPSRural()) ≈ cwind
            @test vertical_dispersion(1.2, Plume, class, CCPSRural()) ≈ vert
        end

        # Puff dispersion
        knowns = [(ClassA, 0.21287234847910436, 0.6879188103871441),
                (ClassB, 0.16556738215041453, 0.6054498545741843),
                (ClassC, 0.11826241582172466, 0.3869880921561032),
                (ClassD, 0.07095744949303479, 0.17041904657983328),
                (ClassE, 0.047304966328689864, 0.11258170198247626),
                (ClassF, 0.023523465599668385, 0.05588182287353654)]
        @testset "Puff dispersion, stability class $class" for (class,cwind,vert) in knowns
            @test crosswind_dispersion(1.2, Puff, class, CCPSUrban()) ≈ cwind
            @test downwind_dispersion(1.2, Puff, class, CCPSUrban()) ≈ cwind
            @test vertical_dispersion(1.2, Puff, class, CCPSUrban()) ≈ vert
        end
    end

    @testset "Windspeed by powerlaw" begin

        u0, z0, p = 3.0, 1.0, 0.108
        a = DryAir(windspeed=u0, windspeed_height=z0, stability=ClassA)
        s = Scenario(Substance(:null,0,0,0,0,0,0,0,0),Release(0,0,0,0,1.0,0,0,0),a)
        @test _windspeed(s,10,CCPSRural()) == _windspeed(a,10,CCPSRural()) == _windspeed(u0,z0,10,ClassA,CCPSRural())

        urban = [(ClassA, 4.237612633868263),
                 (ClassB, 4.237612633868263),
                 (ClassC, 4.754679577383341),
                 (ClassD, 5.334838230116768),
                 (ClassE, 7.53565929452874),
                 (ClassF, 11.943215116604916)]
        @testset "Windspeed, urban terrain, stability class $class" for (class, ans) in urban
            @test  _windspeed(u0,z0,10,class,CCPSUrban()) ≈ ans
        end

        rural = [(ClassA, 3.5246926648185886),
                 (ClassB, 3.5246926648185886),
                 (ClassC, 3.776776235382502),
                 (ClassD, 4.237612633868263),
                 (ClassE, 6.716163415705019),
                 (ClassF, 10.644401677007265)]
        @testset "Windspeed, rural terrain, stability class $class" for (class, ans) in rural
            @test  _windspeed(u0,z0,10,class,CCPSRural()) ≈ ans
        end

    end

end