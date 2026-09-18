@testset "SLAB steady RHS" begin
    idpf,nxtr,vecs,vars,params,dt = GasDispersion.slab._slab_init_hjet(
        2, 1, 11, 50, 61, 0.017031, 2045.90, 239.57, 0.81, 1170000.0,
        4611.80, 603.00, 2976.01, 0.00, 239.57, 107.87, 0.93, 381.0,
        0.00, 1.00, 10.00, 2800.0, [0.0, 1.0, 0.0, 0.0], 0.003, 2.0,
        4.5, 306.2, 21.3, 0.0, 0.0221)

    phase = GasDispersion.slab._slab_steady_phase_state(vecs, vars)
    context = GasDispersion.slab._slab_steady_context(params, phase, idpf, vars)
    projected, entrainment = GasDispersion.slab._slab_steady_project(context.y0, context, vecs.x[1])
    rhs = GasDispersion.slab._slab_steady_rhs(context.y0, context, vecs.x[1])
    next_state, next_entrainment = GasDispersion.slab._slab_steady_ode_step(
        phase, params, idpf, vecs.x[1], 0.001; rmi=vars.rmi, alfg=vars.alfg,
        sru0=vars.sru0, bbx=vars.bbx)
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
    @test all(isfinite, (next_state.r, next_state.h, next_state.rho, next_state.t))
    @test all(isfinite, (next_entrainment.w, next_entrainment.v, next_entrainment.vx))
end