using OrdinaryDiffEqLowOrderRK: RK4

@testset "SLAB steady RHS" begin
    idpf,nxtr,vecs,vars,params,dt = GasDispersion.slab._slab_init_hjet(
        2, 1, 11, 50, 61, 0.017031, 2045.90, 239.57, 0.81, 1170000.0,
        4611.80, 603.00, 2976.01, 0.00, 239.57, 107.87, 0.93, 381.0,
        0.00, 1.00, 10.00, 2800.0, [0.0, 1.0, 0.0, 0.0], 0.003, 2.0,
        4.5, 306.2, 21.3, 0.0, 0.0221)

    phase = GasDispersion.slab._slab_steady_phase_state(vecs, vars)
    context = GasDispersion.slab._slab_steady_context(params, phase, idpf, vars;
        x0=vecs.x[1])
    projected, entrainment = GasDispersion.slab._slab_steady_project(context.y0, context, vecs.x[1])
    rhs = GasDispersion.slab._slab_steady_rhs(context.y0, context, vecs.x[1])
    controls = GasDispersion.slab.SLAB_Steady_Controls(
        vars.rmi, vars.alfg, vars.sru0, vars.bbx)
    input = GasDispersion.slab.SLAB_Steady_IntegratorInput(
        params, phase, idpf, vecs.x[1], vecs.x[1] + 0.001, controls,
        (;dt=0.001, adaptive=false))
    intctx = GasDispersion.slab._slab_steady_integrator(RK4(), input)
    step = GasDispersion.slab.SLAB_Steady_StepInput(input, nothing, nothing)
    next_result = GasDispersion.slab._slab_steady_step!(intctx, step)
    legacy_rhs = zeros(Float64, 11)
    GasDispersion.slab._slab_sub_slope!(legacy_rhs, params, vecs.rho[1], vecs.x[1],
        vecs.h[1], vecs.v[1], vecs.w[1], vecs.b[1], vecs.bb[1], vecs.vg[1],
        vecs.u[1], vecs.wc[1], vecs.cm[1], vars.ft, vars.fu, vars.fv, vars.fw,
        params.othr.bse)

    @test projected.r ≈ phase.r
    @test projected.t ≈ phase.t
    @test projected.rho ≈ phase.rho
    @test rhs ≈ legacy_rhs rtol=1e-5
    @test all(isfinite, rhs)
    @test all(isfinite, (entrainment.w, entrainment.v, entrainment.vx))
    @test all(isfinite, (next_result.state.r, next_result.state.h,
                         next_result.state.rho, next_result.state.t))
    @test isfinite(next_result.cv)

    input = GasDispersion.SLAB_Input(idspl=2,ncalc=1,wms=0.017031,cps=2045.90,
        tbp=239.57,cmed0=0.81,dhe=1170000.0,cpsl=4611.80,rhosl=603.00,
        spb=2976.01,spc=0.00,ts=239.57,qs=107.87,as=0.93,tsd=381.0,
        qtis=0.00,hs=1.00,tav=10.00,xffm=2800.00,zp=[0.0,1.0,0.0,0.0],
        z0=0.003,za=2.0,ua=4.5,ta=306.2,rh=21.3,stab=0.0,ala=0.0221)
    legacy_output = GasDispersion.slab.slab_main(input)
    ode_output = GasDispersion.slab.slab_main(input,RK4();
        steady_solver_kwargs=(;adaptive=false))

    @test legacy_output.steady.ode_solution === nothing
    @test ode_output.steady.ode_solution !== nothing
    @test ode_output.steady.ode_solution.t[1] ≈ vecs.x[1]
    @test ode_output.steady.ode_solution.t[end] > ode_output.steady.ode_solution.t[1]
    @test ode_output.interpolations.ode_solution === ode_output.steady.ode_solution
    @test ode_output.interpolations.cc isa GasDispersion.slab.SLAB_ODE_FieldInterpolation
    @test all(isfinite, (ode_output.interpolations.cc(1.05),
                         ode_output.interpolations.b(1.05),
                         ode_output.interpolations.betac(1.05),
                         ode_output.interpolations.zc(1.05),
                         ode_output.interpolations.sig(1.05)))
end