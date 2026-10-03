function _slab_int_steady_state_impl!(vecs::SLAB_Vecs{F,A},vars::SLAB_Loop_Init{I,F},
                                      params::SLAB_Params{I,F,A},idpf::I,nxtr::I;
                                      solver=SLABLegacySolver(),solver_kwargs=(;)) where {
                                  I <: Integer, F <: AbstractFloat, A <: AbstractVector{F}}
                             
    # unpack parameters
    cmdaa = params.met.cmdaa
    qs = params.spl.qs
    rhoa = params.met.rhoa
    tsd = params.spl.tsd
    qtcs = params.spl.qtcs
    tgon = params.othr.tgon
    bse = params.othr.bse
    urf = params.othr.urf
    cf0 = params.othr.cf0
    rcf = params.othr.rcf
    afa = params.othr.afa

    # unpack intial loop variables
    wss = zero(F)
    ft = vars.ft
    fu = vars.fu
    fv = vars.fv
    fw = vars.fw
    fug = vars.fug
    alfg = vars.alfg
    sru0 = vars.sru0
    htp = vars.htp0
    rmi = vars.rmi
    bx = bx0 = vars.bx
    bbx = bbx0 = vars.bbx
    bvx = bvx0 = vars.bvx0
    bbvx = bbvx0 = vars.bbvx0
    bxs0 = vars.bxs0
    xcc0 = vars.xcc0

    # initialize other loop variables to zero
    xn = timn = xstr = zero(F)
    betax = zero(F)
    steady_vecs = nothing
    transient_vars = nothing
    transient_dt = zero(F)
    transient_started = false

    msfm = vars.msfm
    mnfm = vars.mnfm
    mffm = vars.mffm
    nxi = vars.nxi
    gam = vars.gam
    nssm = params.xtra.nssm

    n = max(1, nxi)
    x = vecs.x[n]
    zc = vecs.zc[n]
    h = vecs.h[n]
    bb = vecs.bb[n]
    b = vecs.b[n]
    cv = vecs.cv[n]
    rho = vecs.rho[n]
    t = vecs.t[n]
    u = vecs.u[n]
    uab = vecs.uab[n]
    cm = vecs.cm[n]
    cmev = vecs.cmev[n]
    cmw = vecs.cmw[n]
    cmwv = vecs.cmwv[n]
    wc = vecs.wc[n]
    vg = vg0 = vecs.vg[n]
    ug = vecs.ug[n]
    w = vecs.w[n]
    v = vecs.v[n]
    vx = vecs.vx[n]
    tim = vecs.tim[n]
    beta = vecs.beta[n]
    qint = vecs.qint[n]
    reference = SLAB_Steady_Reference_State(vars.bbv0,vars.bv0,zc,vars.r0,
        qint,t,cmev,cm,cmw,cmwv,vars.cp0,h,u,uab,b,bb,rho,vg0,wc,htp,beta,
        vars.ubs20)

    xffm = params.fld.xffm
    nstp = nssm*mnfm
    dx = (gam - 1) * (xffm - vecs.x[msfm])/((gam^nstp) - 1)
    work = SLAB_Steady_Workspace(zeros(F,11), zeros(F,11), zeros(F,11),
                                 zeros(F,3), zeros(F,4))
    base = _slab_steady_loop_state(reference.r,reference.bbv,reference.bv,
        reference.zc,reference.qint,h,b,bb,rho,t,u,uab,vg0,vg,wc,htp,w,v,vx,
        cm,cmw,cmwv,cmev,reference.cp,ft,fu,fv,fw,fug,reference.ubs2,beta)
    controls = SLAB_Steady_Controls(rmi, alfg, sru0, bbx)
    integrator_input = SLAB_Steady_IntegratorInput(params, base, idpf, x, x + dx,
        controls, merge((dt=dx,), solver_kwargs))
    integrator = _slab_steady_integrator(solver, integrator_input)
    ode_segments = _slab_steady_ode_segments(integrator, F)

    for nx in nxi:mffm
        for ns in 1:nssm
            xn = x + dx
            base = _slab_steady_loop_state(reference.r,reference.bbv,reference.bv,
                reference.zc,reference.qint,h,b,bb,rho,t,u,uab,vg0,vg,wc,htp,
                w,v,vx,cm,cmw,cmwv,cmev,reference.cp,ft,fu,fv,fw,fug,
                reference.ubs2,beta)
            step_input = SLAB_Steady_IntegratorInput(params, base, idpf, x, xn,
                controls, solver_kwargs)
            step_state = SLAB_Steady_StepState(reference, ug)
            step = SLAB_Steady_StepInput(step_input, step_state, work)
            result = _slab_steady_step!(integrator, step)
            next = result.state
            segment = _slab_steady_ode_segment(integrator,x,xn)
            segment === nothing || push!(ode_segments,segment)
            zc,qint = next.zc,next.qint
            h,b,bb,rho,t,u,uab,vg,wc,htp = next.h,next.b,next.bb,next.rho,next.t,
                next.u,next.uab,next.vg,next.wc,next.htp
            cm,cv,cmw,cmwv,cmev = next.cm,result.cv,next.cmw,next.cmwv,next.cmev
            ft,fu,fv,fw,fug = next.ft,next.fu,next.fv,next.fw,next.fug
            beta,vg0,w,v,vx = next.beta,next.vg0,next.w,next.v,next.vx

            x = xn
            reference, controls = _slab_steady_reference_update(result, controls, rhoa)

            dx = gam*dx

        #660 continue
        end

        _slab_sub_store!(vecs,nx,x,bb,b,vg,cm,t,rho,u,h,cv,beta,w,v,cmdaa,cmw,cmwv,
                         cmev,uab,wc,zc,qint,tim,bbx,bx,betax,ug,vx)
        vecs.tccp[nx] = (qint+qint)/qs

        if qint < 0.5*qtcs
            continue
        else
            idpf = 2
            nxtr = nx
            nxi = nx+1
            dt = dx/u
            steady_vecs = deepcopy(vecs)
            r = 0.25*qs*tsd/cm
            rmi = 0.0
            bbx = r/(rho*bb*h)
            bbx0 = bbx
            vecs.bbx[nx] = bbx
            bbvx = bbx
            bbvx0 = bbx
            bbx = 0.9999*bbx
            bx0 = bbx
            vecs.bx[nx] = bbx
            bvx = bbx
            bvx0 = bbx
            vecs.betax[nx] = sqrt(bbx*bbx-bbx*bbx)/√(3)
            sru0 = r*(u - (1 - cm)*uab)
            fv = bbx*fv
            fu = bbx*fu
            fw = bbx*fw
            fug = bbx*fug
            ft = bbx*ft

            ug = (bb/bbx)*vg
            ug0 = ug
            vecs.ug[nx] = ug

            tvars = SLAB_Loop_Init(nxi,msfm,mnfm,mffm,gam,ft,fu,fv,fw,fug,
                reference.bbv,reference.bv,r,reference.cp,controls.alfg,sru0,
                reference.htp,reference.ubs2,rmi,bx,bbx,bbvx0,bvx0,xcc0,bxs0)

            transient_vars = tvars
            transient_dt = dt
            transient_started = true
            break
        end
    end

    #c   steady state calc of timp
    if nxtr > length(vecs.x)
        xptr = (0.5*params.spl.qtcs/qint)*(vecs.x[end]-vecs.x[1]) + vecs.x[1]
        bxtr = bxs0 + xptr - xcc0
        itr = length(vecs.x)
    else
        xptr = vecs.x[nxtr]
        bxtr = vecs.bbx[nxtr]
        itr = nxtr - 1
    end
    
    if itr > 0
        txt = xptr - xcc0 - xcc0
        txb = txt + xptr
        bxr = (bxtr-bxs0)/(xptr-xcc0)
        for i in 1:itr
            vecs.tim[i] = tsd*(vecs.x[i]+txt)/txb
            vecs.bbx[i] = bxs0 + bxr*(vecs.xccp[i]-xcc0)
            vecs.bx[i] = .9999*vecs.bbx[i]
            vecs.betax[i] = sqrt(vecs.bbx[i]*vecs.bbx[i]-vecs.bx[i]*vecs.bx[i])/√(3)
        end
    end

    # this only runs for an evaporating pool release
    # if idspl == 1
    #     for i in 1:5
    #         vecs.tim[i] = vecs.tim[12-i]
    #     end
    # end
    steady_vecs === nothing && (steady_vecs = deepcopy(vecs))
    transient = transient_started ?
        SLAB_Steady_Transient_Handoff(transient_vars, transient_dt, nxtr) : nothing
    return SLAB_Steady_Phase_Result(steady_vecs, transient,
        _slab_steady_solution(integrator), ode_segments)
end