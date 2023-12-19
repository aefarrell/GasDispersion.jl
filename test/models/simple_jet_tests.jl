@testset "Simple turbulent jet tests" begin
    # example scenario
    Patm = 101325 # Pa
    P1 = 4e5 + Patm # Pa
    T1 = 25 + 273.15 # K
    
    propane = Substance(name = :propane,
                vapor_pressure=0,
                gas_density = 9.7505, # Propane, NIST Webbook
                liquid_density = 526.13, # Propane, NIST Webbook
                reference_temp= T1,
                reference_pressure= P1,
                k=1.15, # heat capacity ratio, from Crane's
                boiling_temp = 231.04, # Propane, NIST Webbook
                latent_heat = 425740, # J/kg, 
                gas_heat_capacity = 1678, # J/kg/K, 
                liquid_heat_capacity = 2520) # J/kg/K
    
    scn = scenario_builder(propane, JetSource; 
           phase = :gas,
           diameter = 0.01, # m
           dischargecoef = 0.85,
           temperature = T1, # K
           pressure = P1,    # Pa
           height = 3.5)     # m, height of hole above the ground
    
    # known answers
    a, b = 6, 5
    ξ² = log(2)/b^2
    x, y, z =  a, √(a^2*ξ²), 1
    c = 0.00014695820366436948

    # horizontal jet
    j = plume(scn, SimpleJet; release_angle=0, k2=a, k3=b)
    @test isa(j, GasDispersion.SimpleJetSolution)
    @test isa(j, Plume)
    @test j(x,y,z) ≈ c
    @test j(-x,y,z) == 0.0

end
