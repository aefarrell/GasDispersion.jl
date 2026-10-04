"""
    _slab_int_steady_state_impl!(vecs, vars, params, idpf, nxtr;
                                 solver=SLABLegacySolver(), solver_kwargs=(;))

Run the steady plume phase using the selected backend. The legacy implementation
retains its fixed spatial stepping schedule; OrdinaryDiffEq integrates the
continuous steady-state system over the full domain and samples output values
without restricting its accepted timesteps to the legacy grid.
"""
function _slab_int_steady_state_impl!(vecs::SLAB_Vecs{F,A},vars::SLAB_Loop_Init{I,F},
                                      params::SLAB_Params{I,F,A},idpf::I,nxtr::I;
                                      solver=SLABLegacySolver(),solver_kwargs=(;)) where {
                                  I <: Integer, F <: AbstractFloat, A <: AbstractVector{F}}
    # Set up the shared initial phase state and backend-neutral controls from
    # the initializer's first available spatial output point.
    qs = params.spl.qs
    tsd = params.spl.tsd
    qtcs = params.spl.qtcs
    cmdaa = params.met.cmdaa
    alfg = vars.alfg
    sru0 = vars.sru0
    rmi = vars.rmi
    bx = bx0 = vars.bx
    bbx = bbx0 = vars.bbx
    bvx = bvx0 = vars.bvx0
    bbvx = bbvx0 = vars.bbvx0
    bxs0 = vars.bxs0
    xcc0 = vars.xcc0
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

    n = max(1, nxi)
    x = vecs.x[n]
    cv = vecs.cv[n]
    ug = vecs.ug[n]
    tim = vecs.tim[n]
    base = _slab_steady_loop_state(_slab_steady_phase_state(vecs, vars, n))
    controls = SLAB_Steady_Controls(rmi, alfg, sru0, bbx)
    integrator_input = SLAB_Steady_IntegratorInput(params, base, idpf, x,
        params.fld.xffm,
        controls, solver_kwargs)
    integrator = _slab_steady_integrator(solver, integrator_input)
    run = _slab_steady_integrate!(integrator, vecs, vars, params, integrator_input;
                                  ug=ug, tim=tim, bbx=bbx, bx=bx, betax=betax)

    # Both backends return the same compact summary; use it to perform the
    # common steady-to-transient handoff and update shared output geometry.
    base, x, dx, cv = run.base, run.x, run.dx, run.base.cv
    nxtr, reference, controls = run.nxtr, run.reference, run.controls
    rmi, alfg, sru0, bbx = controls.rmi, controls.alfg, controls.sru0, controls.bbx

    if run.stopped
        # Preserve the completed steady fields before preparing the live
        # vectors and loop variables for transient continuation.
        steady_vecs = deepcopy(vecs)
        # The transient continuation expands these arrays as it appends each
        # output step, so the adaptive steady segment needs no legacy capacity.
        idpf = 2
        nxi = nxtr+1
        if run.saved_values === nothing
            dt = dx/base.u
        else
            # Choose an initial transient step whose geometric growth spans
            # the remaining domain over the usual number of outer iterations.
            nsteps = max(mffm-vars.nxi+1, 1)
            nssm = params.xtra.nssm
            growth = gam^nssm
            growth_sum = nssm*(growth^nsteps - 1)/(growth - 1)
            dt = (params.fld.xffm-x)/(base.u*growth_sum)
        end
        r = 0.25*qs*tsd/base.cm
        rmi = 0.0
        bbx = r/(base.rho*base.bb*base.h)
        bbx0 = bbx
        vecs.bbx[nxtr] = bbx
        bbvx = bbx
        bbvx0 = bbx
        bbx = 0.9999*bbx
        bx0 = bbx
        vecs.bx[nxtr] = bbx
        bvx = bbx
        bvx0 = bbx
        vecs.betax[nxtr] = sqrt(bbx*bbx-bbx*bbx)/√(3)
        sru0 = r*(base.u - (1 - base.cm)*base.uab)
        fv = bbx*base.fv
        fu = bbx*base.fu
        fw = bbx*base.fw
        fug = bbx*base.fug
        ft = bbx*base.ft

        ug = (base.bb/bbx)*base.vg
        ug0 = ug
        vecs.ug[nxtr] = ug

        tvars = SLAB_Loop_Init(nxi,msfm,mnfm,mffm,gam,ft,fu,fv,fw,fug,
            reference.bbv,reference.bv,r,reference.cp,controls.alfg,sru0,
            reference.htp,reference.ubs2,rmi,bx,bbx,bbvx0,bvx0,xcc0,bxs0)

        transient_vars = tvars
        transient_dt = dt
        transient_started = true
    end

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
        # Backfill release-time and crosswind-width values using the original
        # steady-to-transient geometry convention.
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
    !run.stopped && (steady_vecs = deepcopy(vecs))

    transient = transient_started ?
        SLAB_Steady_Transient_Handoff(transient_vars, transient_dt, nxtr) : nothing
    return SLAB_Steady_Phase_Result(steady_vecs, transient, run.ode_solution,
        run.saved_values)
end
