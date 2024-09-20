module GasDispersionClapeyronExt

using GasDispersion
using Clapeyron

struct ClapeyronSubstance{N<:Union{AbstractString,Symbol},M,F<:Number} <: AbstractSubstance
    name::N
    model::M
    MW::F       # molar weight, kg/mol
    T_ref::F    # reference temperature for densities, K (default 15°C)
    P_ref::F    # reference pressure for densities, Pa (default 1atm)
    k::F        # heat capacity ratio Cp/Cv, unitless (default 1.4)
    T_b::F      # normal boiling point, K
end

function ClapeyronSubstance(m::Clapeyron.EoSModel;reference_temp=288.15,reference_pressure=101325)
    if length(m.components) > 1
        error("GasDispersion does not currently support multicomponent releases :(")
    end

    name = m.components[1] # I assume this is always a Vector{String}
    MW = Clapeyron.mw(m)[1]/1000 # Clapeyron [g/mol] -> GasDispersion [kg/mol]
    T_ref = reference_temp # K
    P_ref = reference_pressure # Pa
    T_b = Clapeyron.saturation_temperature(m,101325)[1] # Clapeyron [K] -> GasDispersion [K]
    k = 1.4

    return ClapeyronSubstance(name,m,promote(MW,T_ref,P_ref,k,T_b)...)
end

GasDispersion.Substance(m::Clapeyron.EoSModel;reference_temp=288.15,reference_pressure=101325) = ClapeyronSubstance(m;reference_temp=reference_temp,reference_pressure=reference_pressure)

Base.isapprox(a::ClapeyronSubstance, b::ClapeyronSubstance) = all([
    getproperty(a,k)≈getproperty(b,k) for k in fieldnames(typeof(a))
    if typeof(getproperty(a,k))<:Number ])

function Base.show(io::IO, mime::MIME"text/plain", s::S) where { S<:ClapeyronSubstance}
    s_type = split(string(S),"{")[1]
    print(io, "$s_type: $(s.name) \n")
    for key in fieldnames(S) 
        if key ∉ [:name,:model]
            val =  getproperty(s, key)
            var =  Symbol(split(string(key),"_")[1])
            unit = GasDispersion.UNITS[var]
            print(io, "    $key: $val $unit \n")
        end
    end
end

# ClapeyronSubstance property getters
GasDispersion._MW(s::ClapeyronSubstance) = s.MW
GasDispersion._vapor_pressure(s::ClapeyronSubstance, T) = Clapeyron.saturation_pressure(s.model, T)[1]

GasDispersion._cp_gas(s::ClapeyronSubstance,T) = 0
GasDispersion._cp_gas(s::ClapeyronSubstance) = GasDispersion._cp_gas(s, s.T_ref)

GasDispersion._cp_liquid(s::ClapeyronSubstance,T) = 0
GasDispersion._cp_liquid(s::ClapeyronSubstance) = GasDispersion._cp_liquid(s, s.T_ref)

GasDispersion._boiling_temperature(s::ClapeyronSubstance) = s.T_b
GasDispersion._latent_heat(s::ClapeyronSubstance,T) = Clapeyron.enthalpy_vap(s.model, T)/s.MW #Clapeyron returns J/mol, expect J/kg
GasDispersion._latent_heat(s::ClapeyronSubstance) = GasDispersion._latent_heat(s, s.T_b)

# density functions
GasDispersion._liquid_density(s::ClapeyronSubstance, T, P) = Clapeyron.mass_density(s.model,P,T; phase=:liquid) # kg/m^3
GasDispersion._liquid_density(s::ClapeyronSubstance) = GasDispersion._liquid_density(s, s.T_ref, s.P_ref)

GasDispersion._gas_density(s::ClapeyronSubstance, T, P) = Clapeyron.mass_density(s.model,P,T; phase=:gas) # kg/m^3
GasDispersion._gas_density(s::ClapeyronSubstance) = GasDispersion._gas_density(s, s.T_ref, s.P_ref)

function GasDispersion._density(s::ClapeyronSubstance, f_l, T, P)
    f_g = 1 - f_l
    ρ_l = _liquid_density(s,T,P)
    ρ_g = _gas_density(s,T,P)

    return 1/(f_l/ρ_l + f_g/ρ_g)
end




end