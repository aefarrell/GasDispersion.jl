# this is probably a bad way of doing things but I want to test some utils
# that are not otherwise exposed to the user
include("../src/utils/utils.jl")

@testset "Monin-Obukhov length tests" begin
    @test monin_obukhov("A", 1.2) ≈ -11.4*1.2^0.10
    @test monin_obukhov("D", 1.2) ≈ Inf
    @test_throws ErrorException monin_obukhov("Q", 1.2)
end

@testset "Pasquill-Gifford dispersion tests" begin
    σy_plume = crosswind_dispersion("A")
    σy_puff = crosswind_dispersion("A"; plume=false)
    @test σy_plume(1.2) ≈ 0.423*((600/600)^0.2)*(1.2^0.9)
    @test σy_puff(1.2) ≈ 0.18*(1.2^0.92)
    @test_throws ErrorException crosswind_dispersion("Q")

    σz_plume = vertical_dispersion("A")
    σz_puff = vertical_dispersion("A"; plume=false)
    @test σz_plume(1.2) ≈ 107.7*(1.2^-1.7172)*exp(0.2770*log(1.2)^2)
    @test σz_puff(1.2) ≈ 0.60*(1.2^0.75)
    @test_throws ErrorException vertical_dispersion("Q")
end

@testset "Windspeed by powerlaw" begin
    u0, z0, p = 3.0, 1.0, 0.108
    u = windspeed(u0, z0, "A")
    @test u(10) ≈ u0*(10/z0)^p
    @test_throws ErrorException windspeed(u0, z0, "Q")
end

@testset "Windspeed by Monin-Obukhov length" begin
    ustar, zR, k = 3.0, 1.0, 0.35
    λ = monin_obukhov("A", zR)
    a(z) = (1-15*(z/λ))^0.25
    Ψ(a) = 2*log((1+a)/2) + log((1+a^2)/2) - 2*atan(a) + π/2

    uA = windspeed(ustar, zR, λ, "A")
    uD = windspeed(ustar, zR, λ, "D")
    uE = windspeed(ustar, zR, λ, "E")

    @test uA(10) ≈ (ustar/k)*(log((10+zR)/zR) - Ψ(a(10)))
    @test uD(10) ≈ (ustar/k)*(log((10+zR)/zR))
    @test uE(10) ≈ (ustar/k)*(log((10+zR)/zR) - 4.7*(10/λ))
    @test_throws ErrorException windspeed(ustar, zR, λ, "Q")
end
