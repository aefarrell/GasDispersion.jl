"""
    _slab_int_steady_state_impl!(vecs, vars, params, idpf, nxtr;
                                 solver=SLABLegacySolver(), solver_kwargs=(;))

Run the steady plume phase using the selected backend. This function owns the
shared SLAB stepping schedule: it prepares each substep's base state and
controls, dispatches a step to the backend, stores the resulting cloud state,
and detects the transition to the transient phase. Backend-specific solution
and interpolation data are collected through dispatch rather than conditionals
in this loop.
"""
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
    alfg = vars.alfg
    sru0 = vars.sru0
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
    cv = vecs.cv[n]
    ug = vecs.ug[n]
    tim = vecs.tim[n]
    base = _slab_steady_loop_state(_slab_steady_phase_state(vecs, vars, n))
    reference = _slab_steady_reference_state(base)

    xffm = params.fld.xffm
    nstp = nssm*mnfm
    dx = (gam - 1) * (xffm - vecs.x[msfm])/((gam^nstp) - 1)
    work = SLAB_Steady_Workspace(zeros(F,11), zeros(F,11), zeros(F,11),
                                 zeros(F,3), zeros(F,4))
    controls = SLAB_Steady_Controls(rmi, alfg, sru0, bbx)
    integrator_input = SLAB_Steady_IntegratorInput(params, base, idpf, x, x + dx,
        controls, merge((dt=dx,), solver_kwargs))
    integrator = _slab_steady_integrator(solver, integrator_input)
    ode_segments = _slab_steady_ode_segments(integrator, F)

    # Each spatial output point can contain several shorter integration steps.
    # Carry the projected phase state forward, while reference and controls
    # track the values that the legacy algorithm uses to initialize the next step.
    for nx in nxi:mffm
        for ns in 1:nssm
            xn = x + dx
            step_input = SLAB_Steady_IntegratorInput(params, base, idpf, x, xn,
                controls, solver_kwargs)
            step_state = SLAB_Steady_StepState(reference, ug)
            step = SLAB_Steady_StepInput(step_input, step_state, work)
            result = _slab_steady_step!(integrator, step)
            next = result.state
            _slab_steady_append_ode_segment!(ode_segments, integrator, x, xn)
            base = _slab_steady_loop_state(next)
            # Cloud volume is derived during the step but is not part of the
            # phase state; retain it separately for storing this output point.
            cv = result.cv
            x = xn
            reference, controls = _slab_steady_reference_update(result, controls, rhoa)

            dx = gam*dx

        #660 continue
        end

        _slab_sub_store!(vecs,nx,x,base.bb,base.b,base.vg,base.cm,base.t,base.rho,
                         base.u,base.h,cv,base.beta,base.w,base.v,cmdaa,base.cmw,
                         base.cmwv,base.cmev,base.uab,base.wc,base.zc,base.qint,
                         tim,bbx,bx,betax,ug,base.vx)
        vecs.tccp[nx] = (base.qint+base.qint)/qs

        if base.qint < 0.5*qtcs
            continue
        else
            # Handoff data captures the steady-to-transient boundary conditions.
            idpf = 2
            nxtr = nx
            nxi = nx+1
            dt = dx/base.u
            steady_vecs = deepcopy(vecs)
            r = 0.25*qs*tsd/base.cm
            rmi = 0.0
            bbx = r/(base.rho*base.bb*base.h)
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
            sru0 = r*(base.u - (1 - base.cm)*base.uab)
            fv = bbx*base.fv
            fu = bbx*base.fu
            fw = bbx*base.fw
            fug = bbx*base.fug
            ft = bbx*base.ft

            ug = (base.bb/bbx)*base.vg
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
        xptr = (0.5*params.spl.qtcs/base.qint)*(vecs.x[end]-vecs.x[1]) + vecs.x[1]
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