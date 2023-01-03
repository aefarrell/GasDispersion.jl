# this is probably a bad way of doing things but I want to test some utils
# that are not otherwise exposed to the user
include("../src/utils/utils.jl")

@testset "Monin-Obukhov length tests" begin
    @test _monin_obukhov(1.2, ClassA) ≈ -11.4*1.2^0.10
    @test _monin_obukhov(1.2, ClassD) ≈ Inf
end

@testset "Pasquill-Gifford dispersion tests" begin
    @test crosswind_dispersion(1.2, Plume, ClassA) ≈ 0.423*((600/600)^0.2)*(1.2^0.9)
    @test crosswind_dispersion(1.2, Puff, ClassA) ≈ 0.18*(1.2^0.92)

    @test vertical_dispersion(1.2, Plume, ClassA) ≈ 107.7*(1.2^-1.7172)*exp(0.2770*log(1.2)^2)
    @test vertical_dispersion(1.2, Puff, ClassA) ≈ 0.60*(1.2^0.75)

end

@testset "Windspeed by powerlaw" begin
    u0, z0, p = 3.0, 1.0, 0.108
    a = Ambient(windspeed=u0, windspeed_height=z0, stability=ClassA)
    s = Scenario(Substance(:null,0,0,0,0,0,0,0,0),Release(0,0,0,0,1.0,0,0,0),a)
    @test _windspeed(a) ≈ u0
    @test _windspeed(a,10) ≈ u0*(10/z0)^p
    @test _windspeed(s) ≈ u0
    @test _windspeed(s,10) ≈ u0*(10/z0)^p
end

@testset "Windspeed by Monin-Obukhov length" begin
    ustar, zR, k = 3.0, 1.0, 0.35
    λ = _monin_obukhov(zR, ClassA)
    a(z) = (1-15*(z/λ))^0.25
    Ψ(a) = 2*log((1+a)/2) + log((1+a^2)/2) - 2*atan(a) + π/2

    @test _windspeed(10, ustar, zR, λ, ClassA) ≈ (ustar/k)*(log((10+zR)/zR) - Ψ(a(10)))
    @test _windspeed(10, ustar, zR, λ, ClassD) ≈ (ustar/k)*(log((10+zR)/zR))
    @test _windspeed(10, ustar, zR, λ, ClassE) ≈ (ustar/k)*(log((10+zR)/zR) - 4.7*(10/λ))
end
