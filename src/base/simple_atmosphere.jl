const DRYAIR = Substance(name = :dryair,
                         molar_weight = 0.028960, # kg/mol
                         vapor_pressure=DIPPRVaporPressure(21.662,-692.39,-0.392,4.7574e-3,1),
                         liquid_density=DIPPRLiquidDensity(0.028960,132.45,2.8963,0.26733,132.45,0.27341),
                         reference_temp=288.15,
                         reference_pressure=101325.0,
                         k=1.4,
                         boiling_temp=78.80,
                         latent_heat=DIPPRLatentHeat(0.028960,132.45,0.8474e7,0.3822,0,0,0),
                         gas_heat_capacity=DIPPRIdealGasHeatCapacity(0.028960,132.45,0.28958e5,0.0939e5,3.012e3,0.0758e5,1484),
                         liquid_heat_capacity=DIPPRLiquidHeatCapacity(Eq1,0.028960,132.45,-214_460,9_185.1,-106.12,0.41616,0.0))

const WATER = Substance(name = :water,
                        molar_weight = 0.018015, # kg/mol
                        vapor_pressure= DIPPRVaporPressure(73.649,-7_258.2,-7.3037,4.1653e-6,2),
                        liquid_density= (T,P) -> 0.018015*1000*(17.863 + 58.606*(T/647.096)^(0.35) - 95.396*(T/647.096)^(2/3) + 213.89*(T/647.096) - 141.26*(T/647.096)^(4/3)),
                        reference_temp=288.15,
                        reference_pressure=101325.0,
                        k=1.3,
                        boiling_temp=373.15,
                        latent_heat=DIPPRLatentHeat(0.018015,647.096,5.2053e7,0.3199,-0.212,0.25795,0.0),
                        gas_heat_capacity=DIPPRIdealGasHeatCapacity(0.018015,647.096,0.33363e5,0.2679e5,2.6105e3,0.08896e5,1169),
                        liquid_heat_capacity=DIPPRLiquidHeatCapacity(Eq1,0.018015,647.096,276_370,-2_090.1,8.125,-0.014116,9.3701e-6))

"""
    SimpleAtmosphere{<:Number,<:StabilityClass}(kwargs...)<:Atmosphere

A simple model of the atmosphere.

# Arguments
- `pressure::Number=101325`: atmospheric pressure, Pa
- `temperature::Number=298.15`: ambient pressure, K
- `windspeed::Number=1.5`: windspeed at anemometer height, m/s
- `windspeed_height::Number=10`: anemometer height, m
- `rel_humidity::Number=0`: relative humidity, %
- `stability::StabilityClass=ClassF()`: Pasquill-Gifford stability class

"""
struct SimpleAtmosphere{F<:Number,S<:StabilityClass} <: Atmosphere
    P::F  # atmospheric pressure, Pa
    T::F  # atmospheric temperature, K
    u::F  # windspeed at windspeed height, m/s
    h::F  # reference height for windspeed, m
    rh::F # relative humidity, %
    stability::S # Pasquill-Gifford stability class
end
SimpleAtmosphere(P,T,u,h,rh,stability) = SimpleAtmosphere(promote(P,T,u,h,rh,)...,stability)
SimpleAtmosphere(; pressure=101325,temperature=298.15,windspeed=1.5,windspeed_height=10,
                   rel_humidity=0.0,stability=ClassF()) = SimpleAtmosphere(pressure,
                   temperature,windspeed,windspeed_height,rel_humidity,stability)

_lapse_rate(a::SimpleAtmosphere{<:Number,<:ClassE}) = 0.020
_lapse_rate(a::SimpleAtmosphere{<:Number,<:ClassF}) = 0.035
_lapse_rate(a::SimpleAtmosphere) = nothing
_rel_humidity(a::SimpleAtmosphere) = a.rh
_surface_roughness(a::SimpleAtmosphere) = 1.0

function _density(a::SimpleAtmosphere, T, P)
    ρ_a = _gas_density(DRYAIR,T,P)
    ρ_w = _gas_density(WATER,T,P)

    rh = _rel_humidity(a)/100
    y_w = rh*min(_vapor_pressure(WATER,T)/P,1)
    y_a = 1-y_w
    return y_a*ρ_a + y_w*ρ_w
end