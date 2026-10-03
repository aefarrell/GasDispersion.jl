"""
Legacy fixed-step RK4 backend for the steady-state phase.

The shared steady loop calls these methods through dispatch. This backend has no
continuous trajectory to retain, so its solution hooks return `nothing`; its
output interpolation is built from the stored SLAB vectors using Akima
interpolation.
"""

"""
Initialize the legacy backend.

The original RK4 algorithm has no persistent integrator state, so initialization
only returns the dispatch token used by the subsequent step calls.
"""
function _slab_steady_integrator(::SLABLegacySolver, input::SLAB_Steady_IntegratorInput)
    return SLABLegacyIntegrator()
end

_slab_steady_append_ode_segment!(segments, ::SLABLegacyIntegrator, x0, x1) = nothing

_slab_steady_ode_segments(::SLABLegacyIntegrator, ::Type{F}) where {F} = nothing
_slab_steady_solution(::SLABLegacyIntegrator) = nothing

"""
Build the legacy output interpolations from tabulated cloud data.

Spatial fields are sorted and interpolated against `x`; cloud-center and
crosswind fields are sorted and interpolated against time.
"""
function _slab_legacy_interpolations(cc::SLAB_CC_Vecs)
    xperm = sortperm(cc.x)
    tperm = sortperm(cc.t)
    return SLAB_Interpolations(
        AkimaInterpolation(cc.cc[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.b[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.betac[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.zc[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.sig[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.xc[tperm], cc.t[tperm]),
        AkimaInterpolation(cc.bx[tperm], cc.t[tperm]),
        AkimaInterpolation(cc.betax[tperm], cc.t[tperm]),
        AkimaInterpolation(cc.bx[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.bbx[xperm], cc.x[xperm]),
        AkimaInterpolation(cc.tcld[xperm], cc.x[xperm]),nothing)
end

_slab_output_interpolations(cc::SLAB_CC_Vecs, ::Nothing, params) =
    _slab_legacy_interpolations(cc)

"""Construct the backwards-compatible output wrapper using legacy interpolations."""
SLAB_Output(params, state, cc) = SLAB_Output(params, state, cc,
    _slab_legacy_interpolations(cc), nothing, nothing)

function _slab_output_interpolations(cc::SLAB_CC_Vecs,
        ::SLAB_Steady_Solution{I,F,A,Nothing,Nothing}, params) where {I,F,A}
    return _slab_legacy_interpolations(cc)
end

"""
Advance one legacy steady-state substep with classical four-stage RK4.

The substep evaluates the legacy slope, solve, thermodynamic, evaluation, and
entrainment routines at each RK stage. Its scratch arrays are reused through
`step.workspace`; the final full phase state and derived cloud volume are
returned together.
"""
function _slab_steady_step!(::SLABLegacyIntegrator, step::SLAB_Steady_StepInput)
    input = step.input
    base, params, idpf = input.base, input.params, input.idpf
    x, xf = input.x0, input.x1
    controls = input.controls
    state = step.state.reference
    work = step.workspace
    f = work.f
    sums = work.sum
    dy = work.dy
    dxxi = work.dxxi
    dxrk = work.dxrk
    fill!(sums, zero(eltype(sums)))
    dx = xf - x
    dxxi .= (dx / 2, dx / 2, dx)
    dxrk .= (dx / 6, dx / 3, dx / 3, dx / 6)

    bbv = bv = qint = zc = r = g = gw = sft = sfu = sfy = sfz = zero(eltype(f))
    cm = base.cm
    cv = cmw = cmwv = cmev = t = _cp = zero(eltype(f))
    rho = base.rho
    u, uab, b, bb, beta, h = base.u, base.uab, base.b, base.bb, base.beta, base.h
    vg, vg0, wc, htp = base.vg, state.vg, base.wc, state.htp
    w, v, vx = base.w, base.v, base.vx
    ubs2, fug = zero(eltype(f)), zero(eltype(f))
    ft, fu, fv, fw = base.ft, base.fu, base.fv, base.fw
    xn = x
    for k in 1:4
        _slab_sub_slope!(f, params, rho, x, h, v, w, b,
                 bb, vg, u, wc, cm, ft, fu, fv, fw,
                 params.othr.bse)
        for j in 1:11
            sums[j] += dxrk[k] * f[j]
            if k == 4
                dy[j] = sums[j]
            else
                dy[j] = dxxi[k] * f[j]
                xn = x + dxxi[k]
            end
        end
        bbv,bv,qint,zc,r,g,gw,sft,sfu,sfy,sfz = _slab_sub_solve(
            params, dy, state.bbv, state.bv, state.zc, state.r, h, state.qint)
        cm,cv,cmw,cmwv,cmev,t,rho,_cp = _slab_sub_thermo(
            params, idpf, xn, zero(eltype(f)), controls.rmi, state.t, state.cmev,
            state.cm, state.cmw, state.cmwv, state.cp, r, state.r, sft,
            params.othr.bse)
        u,uab,b,bb,beta,h,zc,vg,vg0,wc,htp = _slab_sub_eval(
            params, xn, controls.alfg, controls.sru0, zc, state.h, state.u, state.uab,
            state.b, state.bb, r, state.r, bv, state.bv, bbv, state.bbv,
            rho, state.rho, vg0, state.wc, cm, state.htp, b, htp,
            uab, state.beta, vg, wc, h, u, bb, sfu, sfz, sfy, g, gw, params.othr.bse)
        w,v,vx,ubs2,fug,ft,fu,fv,fw = _slab_sub_entran(
            params, idpf, xn, zero(eltype(f)), zero(eltype(f)), zero(eltype(f)),
            state.ubs2, u, step.state.ug, vg, uab, rho, zc, t, h, htp, bb, controls.bbx,
            wc, _cp, params.othr.tgon, params.othr.bse, params.othr.urf,
            params.othr.rcf, params.othr.afa)
    end
    next = SLAB_Steady_Phase_State(r,bbv,bv,g,gw,sft,sfu,sfy,sfz,zc,qint,h,b,bb,rho,t,
        u,uab,beta,vg0,vg,wc,htp,w,v,vx,cm,cmw,cmwv,cmev,_cp,ft,fu,fv,fw,fug,ubs2)
    return SLAB_Steady_StepResult(next, cv)
end
