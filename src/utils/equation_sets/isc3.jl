# Correlations from the ISC3 *User's Guide for the Industrial Source Complex (ISC3) Dispersion Models, Vol II Description of Model Algorithms*

# Correlations for urban terrain
struct ISC3UrbanWind <: PowerLawWind end
ISC3Urban = BasicEquationSet{ISC3UrbanWind,Nothing,BriggsUrbanσy,BriggsUrbanσz}

"""
    windspeed(u0::Number,z0::Number,z::Number, stability, ISC3UrbanWind)
returns the windspeed at height `z` for a given Pasquill-Gifford
stability class, `z` is assumed to be in meters and `u` is in m/s

Assumes urban terrain.

# References
+ Irwin, J. S. 1979. "A Theoretical Variation of the Wind Profile Power Law Exponent as a Function of Surface Roughness, and Stability," *Atmospheric Environment, 13: 191-194.
+ EPA. 1995. *User's Guide for the Industrial Source Complex (ISC3) Dispersion Models, vol 2*. United States Environmental Protection Agency EPA-454/B-95-003b
"""
windspeed(u0::Number,z0::Number,z::Number,::Type{ClassA},::Type{ISC3UrbanWind}) = u0*(z/z0)^0.15
windspeed(u0::Number,z0::Number,z::Number,::Type{ClassB},::Type{ISC3UrbanWind}) = u0*(z/z0)^0.15
windspeed(u0::Number,z0::Number,z::Number,::Type{ClassC},::Type{ISC3UrbanWind}) = u0*(z/z0)^0.20
windspeed(u0::Number,z0::Number,z::Number,::Type{ClassD},::Type{ISC3UrbanWind}) = u0*(z/z0)^0.25
windspeed(u0::Number,z0::Number,z::Number,::Type{ClassE},::Type{ISC3UrbanWind}) = u0*(z/z0)^0.30
windspeed(u0::Number,z0::Number,z::Number,::Type{ClassF},::Type{ISC3UrbanWind}) = u0*(z/z0)^0.30


# Correlations for rural terrain
struct ISC3Ruralσy <: DispersionFunction end
struct ISC3Ruralσz <: DispersionFunction end
ISC3Rural = BasicEquationSet{IrwinRural,Nothing,ISC3Ruralσy,ISC3Ruralσz}


"""
    crosswind_dispersion(x, StabilityClass, ISC3Ruralσy)

Plume crosswind dispersion correlations, for rural terrain

# References
+ EPA. 1995. *User's Guide for the Industrial Source Complex (ISC3) Dispersion Models, vol 2*. United States Environmental Protection Agency EPA-454/B-95-003b
"""
crosswind_dispersion(x::Number, ::Type{ClassA}, ::Type{ISC3Ruralσy}) = 465.11628(x/1000)*tan(0.017453293*(24.1670-2.5334log(x/1000)))
crosswind_dispersion(x::Number, ::Type{ClassB}, ::Type{ISC3Ruralσy}) = 465.11628(x/1000)*tan(0.017453293*(18.3330-1.8096log(x/1000)))
crosswind_dispersion(x::Number, ::Type{ClassC}, ::Type{ISC3Ruralσy}) = 465.11628(x/1000)*tan(0.017453293*(12.5000-1.0857log(x/1000)))
crosswind_dispersion(x::Number, ::Type{ClassD}, ::Type{ISC3Ruralσy}) = 465.11628(x/1000)*tan(0.017453293*(8.3330-0.72382log(x/1000)))
crosswind_dispersion(x::Number, ::Type{ClassE}, ::Type{ISC3Ruralσy}) = 465.11628(x/1000)*tan(0.017453293*(6.2500-0.54287log(x/1000)))
crosswind_dispersion(x::Number, ::Type{ClassF}, ::Type{ISC3Ruralσy}) = 465.11628(x/1000)*tan(0.017453293*(4.1667-0.36191log(x/1000)))


"""
    vertical_dispersion(x, StabilityClass, ISC3Ruralσz)

Plume vertical dispersion correlations, for rural terrain

References:
+ EPA. 1995. *User's Guide for the Industrial Source Complex (ISC3) Dispersion Models, vol 2*. United States Environmental Protection Agency EPA-454/B-95-003b
"""
function vertical_dispersion(x::Number, ::Type{ClassA}, ::Type{ISC3Ruralσz})
    x_km = x/1000 # correlation is in km
    if x_km < 0.10
        a, b = 122.800, 0.94470
        return a*x_km^b
    elseif x_km ≤ 0.15
        a, b = 158.080, 1.05420
        return a*x_km^b
    elseif x_km ≤ 0.20
        a, b = 170.220, 1.09320
        return a*x_km^b
    elseif x_km ≤ 0.25
        a, b = 179.520, 1.12620
        return a*x_km^b
    elseif x_km ≤ 0.30
        a, b = 217.410, 1.26440
        return a*x_km^b
    elseif x_km ≤ 0.40
        a, b = 258.890, 1.40940
        return a*x_km^b
    elseif x_km ≤ 0.50
        a, b = 346.750, 1.72830
        return a*x_km^b
    else
        a, b = 453.850, 2.11660
        sz = a*x_km^b
        return min(sz,5000)
    end
end

function vertical_dispersion(x::Number, ::Type{ClassB}, ::Type{ISC3Ruralσz})
    x_km = x/1000 # correlation is in km
    if x_km ≤ 0.20
        a, b = 90.673, 0.93198
        return a*x_km^b
    elseif x_km ≤ 0.40
        a, b = 98.483, 0.98332
        return a*x_km^b
    else
        a, b = 109.300, 1.09710
        sz = a*x_km^b
        return min(sz,5000)
    end
end

function vertical_dispersion(x::Number, ::Type{ClassC}, ::Type{ISC3Ruralσz})
    x_km = x/1000 # correlation is in km
    sz = 61.141x_km^0.91465
    return min(sz, 5000)
end

function vertical_dispersion(x::Number, ::Type{ClassD}, ::Type{ISC3Ruralσz})
    x_km = x/1000 # correlation is in km
    if x_km ≤ 0.30
        a, b = 34.459, 0.86974
        return a*x_km^b
    elseif x_km ≤ 1.00
        a, b = 32.093, 0.81066
        return a*x_km^b
    elseif x_km ≤ 3.00
        a, b = 32.093, 0.64403
        return a*x_km^b
    elseif x_km ≤ 10.00
        a, b = 33.504, 0.60486
        return a*x_km^b
    elseif x_km ≤ 30.00
        a, b = 36.650, 0.56589
        return a*x_km^b
    else
        a, b = 44.053, 0.51179
        return a*x_km^b
    end
end

function vertical_dispersion(x::Number, ::Type{ClassE}, ::Type{ISC3Ruralσz})
    x_km = x/1000 # correlation is in km
    if x_km ≤ 0.10
        a, b = 24.260, 0.83660
        return a*x_km^b
    elseif x_km ≤ 0.30
        a, b = 23.331, 0.81956
        return a*x_km^b
    elseif x_km ≤ 1.00
        a, b = 21.628, 0.75660
        return a*x_km^b
    elseif x_km ≤ 2.00
        a, b = 21.628, 0.63077
        return a*x_km^b
    elseif x_km ≤ 4.00
        a, b = 22.534, 0.57154
        return a*x_km^b
    elseif x_km ≤ 10.00
        a, b = 24.703, 0.50527
        return a*x_km^b
    elseif x_km ≤ 20.00
        a, b = 26.970, 0.46713
        return a*x_km^b
    elseif x_km ≤ 40.00
        a, b = 35.420, 0.37615
        return a*x_km^b
    else
        a, b = 47.618, 0.29592
        return a*x_km^b
    end
end

function vertical_dispersion(x::Number, ::Type{ClassF}, ::Type{ISC3Ruralσz})
    x_km = x/1000 # correlation is in km
    if x_km ≤ 0.20
        a, b = 15.209, 0.81558
        return a*x_km^b
    elseif x_km ≤ 0.70
        a, b = 14.457, 0.78407
        return a*x_km^b
    elseif x_km ≤ 1.00
        a, b = 13.953, 0.68465
        return a*x_km^b
    elseif x_km ≤ 2.00
        a, b = 13.953, 0.63227
        return a*x_km^b
    elseif x_km ≤ 3.00
        a, b = 14.823, 0.54503
        return a*x_km^b
    elseif x_km ≤ 7.00
        a, b = 16.187, 0.46490
        return a*x_km^b
    elseif x_km ≤ 15.00
        a, b = 17.836, 0.41507
        return a*x_km^b
    elseif x_km ≤ 30.00
        a, b = 22.651, 0.32681
        return a*x_km^b
    elseif x_km ≤ 60.00
        a, b = 27.074, 0.27436
        return a*x_km^b
    else
        a, b = 34.219, 0.21716
        return a*x_km^b
    end
end