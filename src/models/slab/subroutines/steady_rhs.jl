"""Pure steady-state phase state and RHS for the OrdinaryDiffEq backend."""

"""
Complete phase state carried between steady integration substeps.

Besides the 11 values integrated by the ODE solver, this holds the derived
thermodynamic, geometric, and entrainment quantities needed by the next step.
The `g`, `gw`, `sft`, `sfu`, `sfy`, and `sfz` fields are per-step accumulators
and are cleared when preparing a new loop base state.
"""
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

"""
Mutable parameters for evaluating the steady RHS over the current substep.

OrdinaryDiffEq reuses one context while the shared loop changes the interval's
base state and controls. Segment snapshots are made before those fields change
so dense-output interpolation can reproduce the corresponding substep.
"""
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

"""Control values shared by the steady RHS and updated between substeps."""
struct SLAB_Steady_Controls{F <: AbstractFloat}
    rmi::F
    alfg::F
    sru0::F
    bbx::F
end

"""
Backend-neutral inputs for one steady integration interval.

`base` is the current full phase state; `controls` and `idpf` provide the
remaining model configuration, while `x0`/`x1` delimit this interval.
"""
struct SLAB_Steady_IntegratorInput{P,S,I <: Integer,F <: AbstractFloat,C,K}
    params::P
    base::S
    idpf::I
    x0::F
    x1::F
    controls::C
    solver_kwargs::K
end

"""
Reference values carried by the legacy steady loop.

These are distinct from the projected phase state: the legacy equations retain
selected prior-step reference values (including `vg`) while advancing the
current phase state.
"""
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

"""Reference values and the previous crosswind velocity used to start a step."""
struct SLAB_Steady_StepState{R,F <: AbstractFloat}
    reference::R
    ug::F
end

"""Reusable scratch arrays for one legacy RK4 step."""
struct SLAB_Steady_Workspace{F <: AbstractFloat}
    f::Vector{F}
    sum::Vector{F}
    dy::Vector{F}
    dxxi::Vector{F}
    dxrk::Vector{F}
end

"""Inputs grouped for a backend's single steady integration step."""
struct SLAB_Steady_StepInput{I,S,W}
    input::I
    state::S
    workspace::W
end

"""Result of one backend step: full projected phase state and derived cloud volume."""
struct SLAB_Steady_StepResult{S,F <: AbstractFloat}
    state::S
    cv::F
end

function _slab_steady_reference_update(result::SLAB_Steady_StepResult,
                                       controls::SLAB_Steady_Controls, rhoa)
    state = result.state
    reference = _slab_steady_reference_state(state)

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

"""
Outputs of the steady phase before any transient continuation is run.

The shared loop returns the sampled steady vectors, optional transient handoff,
the solver's solution representation, and any native callback storage.
"""
struct SLAB_Steady_Phase_Result{V,T,O,S}
    steady_state::V
    transient::T
    ode_solution::O
    saved_values::S
end

"""
Create the first complete phase state from initialized loop values and vectors.

The reduced ODE accumulators start at zero, while the remaining fields are
initialized from the matching vector index and `SLAB_Loop_Init`.
"""
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

"""
Build the reference state expected by the legacy step equations.

In particular, the prior `vg` field is used rather than `vg0`; preserving this
mapping keeps both integration backends aligned with the legacy loop.
"""
function _slab_steady_reference_state(state::SLAB_Steady_Phase_State)
    # The legacy loop carries state.vg (not state.vg0) into the next reference.
    return SLAB_Steady_Reference_State(state.bbv,state.bv,state.zc,state.r,
        state.qint,state.t,state.cmev,state.cm,state.cmw,state.cmwv,state.cp,state.h,
        state.u,state.uab,state.b,state.bb,state.rho,state.vg,state.wc,state.htp,
        state.beta,state.ubs2)
end

"""
Copy the current RHS context for later interpolation without copying parameters.

The context itself is mutable and will be reused by the integrator. Its phase
state and scalar controls are immutable values, and the parameter bundle is
shared read-only, so copying the context fields is sufficient to preserve an
independent segment snapshot at substantially lower allocation cost than
`deepcopy`.
"""
function _slab_steady_context_snapshot(context::SLAB_Steady_RHS_Context)
    return SLAB_Steady_RHS_Context(context.params, context.base, context.y0,
        context.x0, context.idpf, context.rmi, context.alfg, context.sru0,
        context.bbx, context.tgon, context.bse, context.urf, context.rcf,
        context.afa)
end

"""Build the mutable RHS context and initial 11-component ODE state."""
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

"""
Prepare a phase state for the next steady substep.

All evolving physical fields are retained; only integration accumulators are
reset because each substep starts its solve relative to a fresh base state.
"""
function _slab_steady_loop_state(state::SLAB_Steady_Phase_State{F}) where {F}
    return SLAB_Steady_Phase_State(state.r,state.bbv,state.bv,
        zero(F),zero(F),zero(F),zero(F),zero(F),zero(F),state.zc,state.qint,
        state.h,state.b,state.bb,state.rho,state.t,state.u,state.uab,state.beta,
        state.vg0,state.vg,state.wc,state.htp,state.w,state.v,state.vx,state.cm,
        state.cmw,state.cmwv,state.cmev,state.cp,state.ft,state.fu,state.fv,
        state.fw,state.fug,state.ubs2)
end

"""
Project the 11-component ODE vector into SLAB's complete phase state.

This applies the existing SLAB solve, thermodynamic, evaluation, and entrainment
calculations. Cloud volume and velocity quantities needed by callers are
returned separately as derived values.
"""
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

"""
Evaluate the spatial derivative for the reduced steady-state ODE.

At the initial position, use the base-state slope to preserve the legacy
initialization behavior; later positions use the projected state and current
entrainment values.
"""
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
