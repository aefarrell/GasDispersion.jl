"""Pure steady-state phase state and RHS for the OrdinaryDiffEq backend."""

struct SLAB_Steady_Phase_State{F <: AbstractFloat}
    r::F
    bbv::F
    bv::F
    g::F
    gw::F
    sft::F
    sfu::F
    sfy::F
    sfz::F
    zc::F
    qint::F
    h::F
    b::F
    bb::F
    rho::F
    t::F
    u::F
    uab::F
    beta::F
    vg0::F
    vg::F
    wc::F
    htp::F
    w::F
    v::F
    vx::F
    cm::F
    cmw::F
    cmwv::F
    cmev::F
    cp::F
    ft::F
    fu::F
    fv::F
    fw::F
    fug::F
    ubs2::F
end

mutable struct SLAB_Steady_RHS_Context{I <: Integer, F <: AbstractFloat, A <: AbstractVector{F}}
    params::SLAB_Params{I,F,A}
    base::SLAB_Steady_Phase_State{F}
    y0::SVector{11,F}
    x0::F
    idpf::I
    rmi::F
    alfg::F
    sru0::F
    bbx::F
    tgon::F
    bse::F
    urf::F
    rcf::F
    afa::F
end

struct SLAB_Steady_ODESegment{F <: AbstractFloat,C}
    x0::F
    x1::F
    context::C
end

struct SLAB_Steady_Controls{F <: AbstractFloat}
    rmi::F
    alfg::F
    sru0::F
    bbx::F
end

struct SLAB_Steady_IntegratorInput{P,S,I <: Integer,F <: AbstractFloat,C,K}
    params::P
    base::S
    idpf::I
    x0::F
    x1::F
    controls::C
    solver_kwargs::K
end

struct SLAB_Steady_Reference_State{F <: AbstractFloat}
    bbv::F
    bv::F
    zc::F
    r::F
    qint::F
    t::F
    cmev::F
    cm::F
    cmw::F
    cmwv::F
    cp::F
    h::F
    u::F
    uab::F
    b::F
    bb::F
    rho::F
    vg::F
    wc::F
    htp::F
    beta::F
    ubs2::F
end

struct SLAB_Steady_StepState{R,F <: AbstractFloat}
    reference::R
    ug::F
end

struct SLAB_Steady_Workspace{F <: AbstractFloat}
    f::Vector{F}
    sum::Vector{F}
    dy::Vector{F}
    dxxi::Vector{F}
    dxrk::Vector{F}
end

struct SLAB_Steady_StepInput{I,S,W}
    input::I
    state::S
    workspace::W
end

struct SLAB_Steady_StepResult{S,F <: AbstractFloat}
    state::S
    cv::F
end

function _slab_steady_reference_update(result::SLAB_Steady_StepResult,
                                       controls::SLAB_Steady_Controls, rhoa)
    state = result.state
    # The legacy loop carries next.vg (not next.vg0) into the next reference.
    reference = SLAB_Steady_Reference_State(state.bbv,state.bv,state.zc,state.r,
        state.qint,state.t,state.cmev,state.cm,state.cmw,state.cmwv,state.cp,state.h,
        state.u,state.uab,state.b,state.bb,state.rho,state.vg,state.wc,state.htp,
        state.beta,state.ubs2)

    alfg = controls.alfg
    srug = 0.0
    if !(state.htp > state.h)
        if state.rho > rhoa
            alfg = 0.25
            srug = 0.5 * alfg * grav * (state.rho - rhoa) * state.bb * state.h^2
        else
            alfg = 0.0
        end
    end
    sru0 = state.r * state.u - state.r * (1 - state.cm) * state.uab + srug
    next_controls = SLAB_Steady_Controls(controls.rmi, alfg, sru0, controls.bbx)
    return reference, next_controls
end

struct SLAB_Steady_Transient_Handoff{V,F <: AbstractFloat,I <: Integer}
    vars::V
    dt::F
    index::I
end

struct SLAB_Steady_Phase_Result{V,T,O,G}
    steady_state::V
    transient::T
    ode_solution::O
    ode_segments::G
end

function _slab_steady_phase_state(vecs::SLAB_Vecs{F}, vars::SLAB_Loop_Init{I,F}, index::Integer=1) where {I,F}
    return SLAB_Steady_Phase_State(
        vars.r0, vars.bbv0, vars.bv0, zero(F), zero(F), zero(F), zero(F), zero(F), zero(F),
        vecs.zc[index], vecs.qint[index], vecs.h[index], vecs.b[index], vecs.bb[index],
        vecs.rho[index], vecs.t[index], vecs.u[index], vecs.uab[index], vecs.beta[index],
        vecs.vg[index], vecs.vg[index],
        vecs.wc[index], vars.htp0, vecs.w[index], vecs.v[index], vecs.vx[index],
        vecs.cm[index], vecs.cmw[index], vecs.cmwv[index],
        vecs.cmev[index], vars.cp0, vars.ft, vars.fu, vars.fv, vars.fw, vars.fug, vars.ubs20)
end

function _slab_steady_context(params, state::SLAB_Steady_Phase_State, idpf, vars;
                               rmi=vars.rmi, alfg=vars.alfg, sru0=vars.sru0,
                               bbx=vars.bbx, x0=zero(typeof(state.r)))
    y0 = SVector{11,eltype((state.r,))}(
        state.r, state.bbv, state.bv, state.g, state.sft, state.sfu,
        state.gw, state.zc, state.qint, state.sfy, state.sfz)
    return SLAB_Steady_RHS_Context(params, state, y0, x0, idpf, rmi, alfg,
        sru0, bbx, params.othr.tgon, params.othr.bse, params.othr.urf,
        params.othr.rcf, params.othr.afa)
end

function _slab_steady_loop_state(r0,bbv0,bv0,zc0,qint0,h,b,bb,rho,t,u,uab,vg0,vg,wc,htp,
                                 w,v,vx,cm,cmw,cmwv,cmev,cp0,ft,fu,fv,fw,fug,ubs20,beta)
    F = typeof(r0)
    return SLAB_Steady_Phase_State(r0,bbv0,bv0,zero(F),zero(F),zero(F),zero(F),zero(F),zero(F),
    zc0,qint0,h,b,bb,rho,t,u,uab,beta,vg0,vg,wc,htp,w,v,vx,cm,cmw,cmwv,cmev,cp0,
        ft,fu,fv,fw,fug,ubs20)
end

function _slab_steady_project(u, p::SLAB_Steady_RHS_Context, x)
    base = p.base
    dy = u .- p.y0
    bbv,bv,qint,zc,r,g,gw,sft,sfu,sfy,sfz = _slab_sub_solve(
        p.params, dy, base.bbv, base.bv, base.zc, base.r, base.h, base.qint)
    cm,cv,cmw,cmwv,cmev,t,rho,cp = _slab_sub_thermo(
        p.params, p.idpf, x, zero(eltype(u)), p.rmi, base.t, base.cmev,
        base.cm, base.cmw, base.cmwv, base.cp, r, base.r, sft, p.bse)
    uvel,uab,b,bb,beta,h,zc,vg,vg0,wc,htp = _slab_sub_eval(
        p.params, x, p.alfg, p.sru0, zc, base.h, base.u, base.uab, base.b,
        base.bb, r, base.r, bv, base.bv, bbv, base.bbv, rho, base.rho,
        base.vg0, base.wc, cm, base.htp, base.b, base.htp, base.uab, base.beta,
        base.vg, base.wc, base.h, base.u, base.bb, sfu, sfz, sfy, g, gw, p.bse)
    w,v,vx,ubs2,fug,ft,fu,fv,fw = _slab_sub_entran(
        p.params, p.idpf, x, zero(eltype(u)), zero(eltype(u)), zero(eltype(u)),
        base.ubs2, uvel, zero(eltype(u)), vg, uab, rho, zc, t, h, htp, bb,
        p.bbx, wc, cp, p.tgon, p.bse, p.urf, p.rcf, p.afa)
    return SLAB_Steady_Phase_State(r,bbv,bv,g,gw,sft,sfu,sfy,sfz,zc,qint,h,b,bb,rho,t,
        uvel,uab,beta,vg0,vg,wc,htp,w,v,vx,cm,cmw,cmwv,cmev,cp,ft,fu,fv,fw,fug,ubs2),
        (cv=cv, vg0=vg0, w=w, v=v, vx=vx)
end

function _slab_steady_rhs(u, p::SLAB_Steady_RHS_Context, x)
    base = p.base
    state, entrainment = _slab_steady_project(u, p, x)
    f = zeros(eltype(u), 11)
    if x == p.x0
        _slab_sub_slope!(f, p.params, base.rho, x, base.h, base.v, base.w,
            base.b, base.bb, base.vg, base.u, base.wc, base.cm, base.ft,
            base.fu, base.fv, base.fw, p.bse)
        return SVector{11}(f)
    end
    _slab_sub_slope!(f, p.params, state.rho, x, state.h, entrainment.v,
        entrainment.w, state.b, state.bb, state.vg, state.u, state.wc, state.cm,
        state.ft, state.fu, state.fv, state.fw, p.bse)
    return SVector{11}(f)
end
