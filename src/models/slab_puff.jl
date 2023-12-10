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

function puff(scenario::Scenario, ::Type{SLAB}, eqs::EquationSet=DefaultSet(); 
              t_av=10, x_max=2000)
    Pᵣ = scenario.substance.P_ref
    Tᵣ = scenario.substance.T_ref
    ρⱼ = _gas_density(scenario.substance)
    R = 8.31446261815324
    MW = ρⱼ*R*Tᵣ/Pᵣ
    c_max = 1.0

    stab = _slab_stab( _stability(scenario) )

    inp = SLAB_Input(;idspl = 2,
                     ncalc = 1,
                     wms = MW,
                     cps = _cp_gas(scenario.substance),
                     tbp = _boiling_temperature(scenario.substance),
                     cmed0 = _release_liquid_fraction(scenario),
                     dhe = _latent_heat(scenario.substance),
                     cpsl = _cp_liquid(scenario.substance),
                     rhosl = _liquid_density(scenario.substance),
                     spb = -1.0,
                     spc = 0.0,
                     ts = _release_temperature(scenario),
                     qs = _mass_rate(scenario),
                     as = _release_area(scenario),
                     tsd = _duration(scenario),
                     qtis = 0.00,
                     hs = _release_height(scenario),
                     tav = t_av,
                     xffm = x_max,
                     zp = [0.0],
                     z0 = 1.0,
                     za = _windspeed_height(scenario),
                     ua = _windspeed(scenario),
                     ta = _atmosphere_temperature(scenario),
                     rh = _rel_humidity(scenario.atmosphere),
                     stab = stab,
                     ala = 0.0)
    out = slab_main(inp)
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