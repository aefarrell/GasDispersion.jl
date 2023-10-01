@testset "ISC3 Equation Sets" begin

    @testset "Pasquill-Gifford dispersion tests" begin

        # Plume dispersion
        urban = [(ClassA, 0.38390787316433456, 0.2881727481910807),
                 (ClassB, 0.38390787316433456, 0.2881727481910807),
                 (ClassC, 0.26393666280048, 0.24),
                 (ClassD, 0.19195393658216728, 0.16796976816235135),
                 (ClassE, 0.13196833140024, 0.09591371646531512),
                 (ClassF, 0.13196833140024, 0.09591371646531512)]
        @testset "Plume dispersion, urban terrain, stability class $class" for (class,cwind,vert) in urban
            @test crosswind_dispersion(1.2, Plume, class, ISC3Urban()) ≈ cwind
            @test vertical_dispersion(1.2, Plume, class, ISC3Urban()) ≈ vert
        end

        rural = [(ClassA, 0.48870396813625905, [(50, 7.246283645973222), (125, 17.6538512508938), (175, 25.32210358392127), (225, 33.46114450376929), (275, 42.49832115580982), (350, 58.955561122372494), (450, 87.22955507375895), (550, 128.0454080208172)]),
                 (ClassB, 0.3288132138806368, [(100, 10.604690180980183), (300, 30.144226325216724), (500, 51.092852947678885)]),
                 (ClassC, 0.20096244824229573, [(100, 7.44187785547111)]),
                 (ClassD, 0.1309207626723034, [(100, 4.651174892531855), (500, 18.29689264165363), (2000, 50.15135417398994), (5000, 88.69020460936578), (15000, 169.67281358959912), (35000, 271.77785930262854)]),
                 (ClassE, 0.09742134712173195, [(50, 1.979015073784176), (200, 6.23857638464594), (500, 12.80138815568348), (1500, 27.931190340632067), (3000, 42.22135548587303), (7500, 68.37414410048738), (15000, 95.55830909365176), (30000, 127.31152395989899), (50000, 151.5410717370677)]),
                 (ClassF, 0.0645858766834503, [(100, 2.3255231110829815), (500, 8.395558503802999), (800, 11.976175562087187), (1500, 18.030377292486545), (2500, 24.42448141869856), (5000, 34.207199596890135), (10000, 46.38392157098271), (20000, 60.29440210746138), (45000, 76.93568233551231), (75000, 87.38876819912205)])]
        @testset "Plume dispersion, rural terrain, stability class $class" for (class,cwind,verts) in rural
            @test crosswind_dispersion(1.2, Plume, class, ISC3Rural()) ≈ cwind
            for (x,sz) in verts
                @test vertical_dispersion(x, Plume, class, ISC3Rural()) ≈ sz
            end
        end

        # Puff dispersion
        # just default behaviour
        @test crosswind_dispersion(1.2, Puff, ClassA, ISC3Rural()) == crosswind_dispersion(1.2, Puff, ClassA, DefaultSet())
        @test vertical_dispersion(1.2, Puff, ClassA, ISC3Rural()) == vertical_dispersion(1.2, Puff, ClassA, DefaultSet())
        @test downwind_dispersion(1.2, Puff, ClassA, ISC3Rural()) == downwind_dispersion(1.2, Puff, ClassA, DefaultSet())

    end

    @testset "Windspeed by powerlaw" begin

        u0, z0 = 3.0, 1.0
        a = DryAir(windspeed=u0, windspeed_height=z0, stability=ClassA)
        s = Scenario(Substance(:null,0,0,0,0,0,0,0,0),Release(0,0,0,0,1.0,0,0,0),a)
        @test _windspeed(s,10,ISC3Rural()) == _windspeed(a,10,ISC3Rural()) == _windspeed(u0,z0,10,ClassA,ISC3Rural())

        urban = [(ClassA, 4.237612633868263),
                 (ClassB, 4.237612633868263),
                 (ClassC, 4.754679577383341),
                 (ClassD, 5.334838230116768),
                 (ClassE, 5.985786944906638),
                 (ClassF, 5.985786944906638)]
        @testset "Windspeed, urban terrain, stability class $class" for (class, ans) in urban
            @test  _windspeed(u0,z0,10,class,ISC3Urban()) ≈ ans
        end

        rural = [(ClassA, 3.5246926648185886),
                 (ClassB, 3.5246926648185886),
                 (ClassC, 3.776776235382502),
                 (ClassD, 4.237612633868263),
                 (ClassE, 6.716163415705019),
                 (ClassF, 10.644401677007265)]
        @testset "Windspeed, rural terrain, stability class $class" for (class, ans) in rural
            @test  _windspeed(u0,z0,10,class,ISC3Rural()) ≈ ans
        end

    end

end