include("slab/slab.jl")

using .slab

# defining type for dispatch
struct SLAB <: PuffModel end

struct SLABSolution{I <: Integer, F <: Number, V <: AbstractVector{F}, S} <: Puff
    scenario::Scenario
    model::Symbol
    in::SLAB_Input{I,F,V}
    out::SLAB_Output{I,F,V}
    c₀::F
    cc::S
    b::S
    betac::S
    zc::S
    sig::S
    xc::S
    bx::S
    betax::S
end

# SLAB stability mapping
_slab_stab(::Type{ClassA}) = 1.0
_slab_stab(::Type{ClassB}) = 2.0
_slab_stab(::Type{ClassC}) = 3.0
_slab_stab(::Type{ClassD}) = 4.0
_slab_stab(::Type{ClassE}) = 5.0
_slab_stab(::Type{ClassF}) = 6.0

# SLAB ala mapping
_slab_ala(::SimpleAtmosphere) = 0.0

# Recover the constants of the Antoine equation
function _slab_antoine(s::Scenario)
    if s.substance.P_v isa Number
        return Antoine(0,-1.0,0)
    elseif s.substance.P_v isa Antoine
        return s.substance.P_v
    else
        # Find three temperatures that cover the range of the scenario
        T1 = _boiling_temperature(s.substance)
        T3 = max(_release_temperature(s),_atmosphere_temperature(s))
        T2 = (T1+T3)/2

        # Curve fit to these three points
        P_v = (T) -> _vapor_pressure(s.substance, T)
        T = [T1,T2,T3]
        P = P_v.(T)

        # Fit the linear form y = a1 + a2*x + a3*x*y
        y = log.(P)
        x = T.^-1
        X = [ ones(typeof(T1),3) x x.*y ]
        a = X\y

        # Recover Antoine coefficients
        A = a[1]
        C = -1*a[3]
        B = A*C - a[2]
        return Antoine(A,B,C)
    end
end

@doc doc"""
    puff(::Scenario, SLAB; kwargs...)

Returns the solution to the SLAB horizontal jet dispersion model for the given
scenario.

# References
+ Ermak, Donald L. 1990. *User's Manual for SLAB: An Atmospheric Dispersion Model For Denser-Than-Air Releases* Lawrence Livermore National Laboratory

# Arguments
- `t_av::Number=10`: averaging time, seconds
- `x_max::Number=2000`: maximum downwind distance, meters, this defines the problem domain

"""
function puff(scenario::Scenario, ::Type{SLAB}, eqs::EquationSet=DefaultSet(); 
              t_av=10, x_max=2000)
    c_max = 1.0
    stab = _slab_stab( _stability(scenario) )
    antoine = _slab_antoine(scenario)
    inp = SLAB_Input(;idspl = 2, # horizontal jet
                     ncalc = 1,
                     wms = _MW(scenario.substance),
                     cps = _cp_gas(scenario.substance),
                     tbp = _boiling_temperature(scenario.substance),
                     cmed0 = _release_liquid_fraction(scenario),
                     dhe = _latent_heat(scenario.substance),
                     cpsl = _cp_liquid(scenario.substance),
                     rhosl = _liquid_density(scenario.substance),
                     spb = antoine.B,
                     spc = antoine.C,
                     ts = _release_temperature(scenario),
                     qs = _mass_rate(scenario),
                     as = _release_area(scenario),
                     tsd = _duration(scenario),
                     qtis = 0.00,
                     hs = _release_height(scenario),
                     tav = t_av,
                     xffm = x_max,
                     zp = [0.0],
                     z0 = _surface_roughness(scenario.atmosphere),
                     za = _windspeed_height(scenario),
                     ua = _windspeed(scenario),
                     ta = _atmosphere_temperature(scenario),
                     rh = _rel_humidity(scenario.atmosphere),
                     stab = stab,
                     ala = _slab_ala(scenario.atmosphere))

    # run SLAB and collect output
    out = slab_main(inp)

    # SLAB can return results out of order
    # since v8.0.1 this causes DataInterpolations to error
    xperm = sortperm(out.cc.x)
    tperm = sortperm(out.cc.t)
    return SLABSolution(scenario,:SLAB,inp,out,c_max,
                        AkimaInterpolation(out.cc.cc[xperm], out.cc.x[xperm]),
                        AkimaInterpolation(out.cc.b[xperm], out.cc.x[xperm]),
                        AkimaInterpolation(out.cc.betac[xperm], out.cc.x[xperm]),
                        AkimaInterpolation(out.cc.zc[xperm], out.cc.x[xperm]),
                        AkimaInterpolation(out.cc.sig[xperm], out.cc.x[xperm]),
                        AkimaInterpolation(out.cc.xc[tperm], out.cc.t[tperm]),
                        AkimaInterpolation(out.cc.bx[tperm], out.cc.t[tperm]),
                        AkimaInterpolation(out.cc.betax[tperm], out.cc.t[tperm]))
end

function puff(scenario::Scenario{<:AbstractSubstance,<:VerticalJet,<:Atmosphere}, ::Type{SLAB}, eqs::EquationSet=DefaultSet(); release_angle::Number=π/2, k2::Number=6.0, k3::Number=5.0)
    c_max = 1.0
    stab = _slab_stab( _stability(scenario) )
    antoine = _slab_antoine(scenario)
    inp = SLAB_Input(;idspl = 3, # vertical jet
                     ncalc = 1,
                     wms = _MW(scenario.substance),
                     cps = _cp_gas(scenario.substance),
                     tbp = _boiling_temperature(scenario.substance),
                     cmed0 = _release_liquid_fraction(scenario),
                     dhe = _latent_heat(scenario.substance),
                     cpsl = _cp_liquid(scenario.substance),
                     rhosl = _liquid_density(scenario.substance),
                     spb = antoine.B,
                     spc = antoine.C,
                     ts = _release_temperature(scenario),
                     qs = _mass_rate(scenario),
                     as = _release_area(scenario),
                     tsd = _duration(scenario),
                     qtis = 0.00,
                     hs = _release_height(scenario),
                     tav = t_av,
                     xffm = x_max,
                     zp = [0.0],
                     z0 = _surface_roughness(scenario.atmosphere),
                     za = _windspeed_height(scenario),
                     ua = _windspeed(scenario),
                     ta = _atmosphere_temperature(scenario),
                     rh = _rel_humidity(scenario.atmosphere),
                     stab = stab,
                     ala = _slab_ala(scenario.atmosphere))
    # run SLAB and collect output
    out = slab_main(inp)

    # SLAB can return results out of order
    # since v8.0.1 this causes DataInterpolations to error
    xperm = sortperm(out.cc.x)
    tperm = sortperm(out.cc.t)
    return SLABSolution(scenario,:SLAB,inp,out,c_max,
                        AkimaInterpolation(out.cc.cc, out.cc.x),
                        AkimaInterpolation(out.cc.b, out.cc.x),
                        AkimaInterpolation(out.cc.betac, out.cc.x),
                        AkimaInterpolation(out.cc.zc, out.cc.x),
                        AkimaInterpolation(out.cc.sig, out.cc.x),
                        AkimaInterpolation(out.cc.xc, out.cc.t),
                        AkimaInterpolation(out.cc.bx, out.cc.t),
                        AkimaInterpolation(out.cc.betax, out.cc.t))

end

function (s::SLABSolution)(x,y,z,t)
    h = s.in.hs
    x_max = s.in.xffm
    c_max = s.c₀
    # domain check
    if (x==0)&&(y==0)&&(z==h)
        return c_max
    elseif (x<0)||(z<0)||(t<0)
        return 0.0
    elseif x ≥ x_max
        error("Outside the domain of the solution, x ≥ $x_max")
    else
        x = x+1.0 # slab assumes the emission point is x=1
        cc = s.cc(x)
        xa = (x - s.xc(t) + s.bx(t))/(√(2)*s.betax(t))
        xb = (x - s.xc(t) - s.bx(t))/(√(2)*s.betax(t))
        ya = (y + s.b(x))/(√(2)*s.betac(x))
        yb = (y - s.b(x))/(√(2)*s.betac(x))
        za = (z - s.zc(x))/(√(2)*s.sig(x))
        zb = (z + s.zc(x))/(√(2)*s.sig(x))

        c = cc*erf(xb,xa)*erf(yb,ya)*(exp(-za^2) + exp(-zb^2))
        return min(c,c_max)
    end
end