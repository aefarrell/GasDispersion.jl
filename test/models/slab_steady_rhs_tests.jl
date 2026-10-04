using OrdinaryDiffEqTsit5: Tsit5
using DataInterpolations: AkimaInterpolation

@testset "SLAB steady RHS" begin
    idpf,nxtr,vecs,vars,params,dt = GasDispersion.slab._slab_init_hjet(
        2, 1, 11, 50, 61, 0.017031, 2045.90, 239.57, 0.81, 1170000.0,
        4611.80, 603.00, 2976.01, 0.00, 239.57, 107.87, 0.93, 381.0,
        0.00, 1.00, 10.00, 2800.0, [0.0, 1.0, 0.0, 0.0], 0.003, 2.0,
        4.5, 306.2, 21.3, 0.0, 0.0221)

    phase = GasDispersion.slab._slab_steady_phase_state(vecs, vars)
    context = GasDispersion.slab._slab_steady_context(params, phase, idpf, vars;
        x0=vecs.x[1])
    # Segment snapshots must preserve mutable context values without duplicating
    # the shared parameter object or its input arrays.
    context_snapshot = GasDispersion.slab._slab_steady_context_snapshot(context)
    context.x0 += 1
    @test context_snapshot !== context
    @test context_snapshot.x0 == vecs.x[1]
    @test context_snapshot.params === context.params
    @test context_snapshot.params.fld.zp === context.params.fld.zp
    context.x0 = vecs.x[1]
    projected, entrainment = GasDispersion.slab._slab_steady_project(context.y0, context, vecs.x[1])
    rhs = GasDispersion.slab._slab_steady_rhs(context.y0, context, vecs.x[1])
    controls = GasDispersion.slab.SLAB_Steady_Controls(
        vars.rmi, vars.alfg, vars.sru0, vars.bbx)
    input = GasDispersion.slab.SLAB_Steady_IntegratorInput(
        params, phase, idpf, vecs.x[1], vecs.x[1] + 0.001, controls,
        (;dt=0.001, adaptive=false))
    intctx = GasDispersion.slab._slab_steady_integrator(
        GasDispersion.slab.SLABLegacySolver(), input)
    workspace = GasDispersion.slab.SLAB_Steady_Workspace(
        zeros(11), zeros(11), zeros(11), zeros(3), zeros(4))
    step = GasDispersion.slab.SLAB_Steady_StepInput(input,
        GasDispersion.slab.SLAB_Steady_StepState(
            GasDispersion.slab._slab_steady_reference_state(phase), vecs.ug[1]),
        workspace)
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
    @test all(isfinite, (next_result.r, next_result.h,
                         next_result.rho, next_result.t))
    @test isfinite(next_result.cv)

    reference, next_controls = GasDispersion.slab._slab_steady_reference_update(
        next_result, controls, params.met.rhoa)
    expected_alfg = next_result.htp > next_result.h ? controls.alfg :
        next_result.rho > params.met.rhoa ? 0.25 : 0.0
    srug = next_result.htp > next_result.h ||
        next_result.rho <= params.met.rhoa ? 0.0 :
        0.5 * expected_alfg * params.xtra.grav * (next_result.rho - params.met.rhoa) *
            next_result.bb * next_result.h^2
    @test reference.vg == next_result.vg
    @test next_controls.alfg == expected_alfg
    @test next_controls.sru0 ≈ next_result.r * next_result.u -
        next_result.r * (1 - next_result.cm) * next_result.uab + srug
    @test next_controls.rmi == controls.rmi
    @test next_controls.bbx == controls.bbx

    input = GasDispersion.SLAB_Input(idspl=2,ncalc=1,wms=0.017031,cps=2045.90,
        tbp=239.57,cmed0=0.81,dhe=1170000.0,cpsl=4611.80,rhosl=603.00,
        spb=2976.01,spc=0.00,ts=239.57,qs=107.87,as=0.93,tsd=381.0,
        qtis=0.00,hs=1.00,tav=10.00,xffm=2800.00,zp=[0.0,1.0,0.0,0.0],
        z0=0.003,za=2.0,ua=4.5,ta=306.2,rh=21.3,stab=0.0,ala=0.0221)
    legacy_output = GasDispersion.slab.slab_main(input)
    ode_output = GasDispersion.slab.slab_main(input,Tsit5())

    # The legacy backend has no dense solution; the ODE backend retains one
    # whose accepted steps are independent of the stored output grid.
    @test legacy_output.steady.ode_solution === nothing
    @test ode_output.steady.ode_solution !== nothing
    @test ode_output.steady isa GasDispersion.slab.SLAB_ODE_Steady_Solution
    @test ode_output.steady.state.x == ode_output.steady.saved_values.t
    @test length(ode_output.steady.state.x) ==
        length(ode_output.steady.saved_values.saveval)
    @test length(ode_output.steady.state.x) !=
        length(legacy_output.steady.state.x)
    @test length(ode_output.s.x) == length(ode_output.cc.x)
    @test ode_output.s.x[end] >= ode_output.p.fld.xffm
    @test ode_output.steady.ode_solution.t[1] ≈ vecs.x[1]
    @test ode_output.steady.ode_solution.t[end] > ode_output.steady.ode_solution.t[1]
    legacy_grid = legacy_output.steady.state.x
    accepted_times = ode_output.steady.ode_solution.t[2:end-1]
    @test any(t -> minimum(abs.(legacy_grid .- t)) >
        1e-8*max(abs(t), 1.0), accepted_times)
    @test ode_output.interpolations.ode_solution === ode_output.steady.ode_solution
    @test ode_output.interpolations.cc isa AkimaInterpolation
    @test ode_output.interpolations.b isa AkimaInterpolation
    @test ode_output.interpolations.betac isa AkimaInterpolation
    @test ode_output.interpolations.sig isa AkimaInterpolation
    legacy_interpolations = GasDispersion.slab._slab_legacy_interpolations(ode_output.cc)
    @test ode_output.interpolations.zc(2500.0) ≈ legacy_interpolations.zc(2500.0)
    saved_values = ode_output.steady.saved_values
    sample = cld(length(saved_values.t), 2)
    sample_x = saved_values.t[sample]
    sample_state = saved_values.saveval[sample]
    sample_cv = sample_state.cv
    sample_bx = legacy_interpolations.bx_x(sample_x)
    sample_bbx = legacy_interpolations.bbx_x(sample_x)
    sample_tcld = legacy_interpolations.tcld(sample_x)
    sample_point = GasDispersion.slab._slab_editcc_point(ode_output.p,sample_x,
        saved_values.t[1],sample_state.zc,sample_state.h,sample_state.b,
        sample_state.beta,sample_state.uab,sample_state.cm,sample_cv,
        sample_bx,sample_bbx,sample_tcld)
    @test ode_output.interpolations.cc(sample_x) ≈ sample_point.cc
    @test ode_output.interpolations.b(sample_x) ≈ sample_state.b
    @test ode_output.interpolations.betac(sample_x) ≈ sample_point.betac
    @test ode_output.interpolations.zc(sample_x) ≈ sample_state.zc
    @test ode_output.interpolations.sig(sample_x) ≈ sample_point.sig
    @test all(isfinite, (ode_output.interpolations.cc(1.05),
                         ode_output.interpolations.b(1.05),
                         ode_output.interpolations.betac(1.05),
                         ode_output.interpolations.zc(1.05),
                         ode_output.interpolations.sig(1.05)))

    handoff = findfirst(qint -> qint >= 0.5*legacy_output.p.spl.qtcs,
        legacy_output.steady.state.qint)
    @test handoff !== nothing
    # Compare both solutions at shared positions, allowing ordinary solver
    # error without requiring the adaptive ODE trajectory to match legacy RK4.
    legacy_x = legacy_output.steady.state.x[2:handoff]
    comparison_indices = filter(i -> legacy_x[1] <=
        ode_output.steady.state.x[i] <= min(legacy_x[end],
            ode_output.steady.ode_solution.t[end]),
        2:length(ode_output.steady.state.x))
    for field in (:rho, :t, :u, :cm, :qint)
        legacy_values = getproperty(legacy_output.steady.state, field)[2:handoff]
        legacy_interpolation = AkimaInterpolation(legacy_values, legacy_x)
        ode_x = ode_output.steady.state.x[comparison_indices]
        ode_values = getproperty(ode_output.steady.state, field)[comparison_indices]
        legacy_values = legacy_interpolation.(ode_x)
        relative_error = maximum(abs.(legacy_values .- ode_values)) /
            maximum(abs.(legacy_values))
        @test relative_error <= 0.2
    end
end