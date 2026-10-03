"""Legacy fixed-step RK4 implementation for the steady-state phase."""

function _slab_steady_integrator(::SLABLegacySolver, base, params, idpf, x, xf;
                                  solver_kwargs=(;), kwargs...)
    return SLABLegacyIntegrator()
end

_slab_steady_ode_segment(::SLABLegacyIntegrator, x0, x1) = nothing

function _slab_steady_step!(::SLABLegacyIntegrator, base, params, idpf, x, xf;
                            work, state, kwargs...)
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
    vg, vg0, wc, htp = base.vg, state.vg0, base.wc, state.htp0
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
            params, dy, state.bbv0, state.bv0, state.zc0, state.r0, h, state.qint0)
        cm,cv,cmw,cmwv,cmev,t,rho,_cp = _slab_sub_thermo(
            params, idpf, xn, zero(eltype(f)), state.rmi, state.t0, state.cmev0,
            state.cm0, state.cmw0, state.cmwv0, state.cp0, r, state.r0, sft,
            params.othr.bse)
        u,uab,b,bb,beta,h,zc,vg,vg0,wc,htp = _slab_sub_eval(
            params, xn, state.alfg, state.sru0, zc, state.h0, state.u0, state.uab0,
            state.b0, state.bb0, r, state.r0, bv, state.bv0, bbv, state.bbv0,
            rho, state.rho0, vg0, state.wc0, cm, state.htp0, b, htp,
            uab, state.beta, vg, wc, h, u, bb, sfu, sfz, sfy, g, gw, params.othr.bse)
        w,v,vx,ubs2,fug,ft,fu,fv,fw = _slab_sub_entran(
            params, idpf, xn, zero(eltype(f)), zero(eltype(f)), zero(eltype(f)),
            state.ubs20, u, state.ug, vg, uab, rho, zc, t, h, htp, bb, state.bbx,
            wc, _cp, params.othr.tgon, params.othr.bse, params.othr.urf,
            params.othr.rcf, params.othr.afa)
    end
    next = SLAB_Steady_Phase_State(r,bbv,bv,g,gw,sft,sfu,sfy,sfz,zc,qint,h,b,bb,rho,t,
        u,uab,beta,vg0,vg,wc,htp,w,v,vx,cm,cmw,cmwv,cmev,_cp,ft,fu,fv,fw,fug,ubs2)
    return next, (cv=cv, vg0=vg0, w=w, v=v, vx=vx)
end
