__precompile__()

module slab

using OrdinaryDiffEq: ODEProblem, solve, init, step!, reinit!, set_proposed_dt!
using OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm
using StaticArrays
using DataInterpolations: AkimaInterpolation

export SLAB_Input, SLAB_Output
export SLABLegacySolver
export slab_main

# defining structs, how the data is passed into and out of SLAB
include("structs.jl")

# declare global variables and constants
include("globals.jl")


# code is organized differently than in SLAB, the functions and subroutines
# are defined first, the main program is the function _slab_main
include("functions.jl")
include("subroutines.jl")

"""
    slab_main(inp::SLAB_Input)

to do doc string
"""
slab_main(inp::SLAB_Input; solver=SLABLegacySolver(), kwargs...) = slab_main(inp, solver; kwargs...)

function slab_main(inp::SLAB_Input, solver; steady_solver_kwargs=(;), kwargs...)
    return slab_main(inp.idspl,inp.ncalc,inp.wms,inp.cps,inp.tbp,inp.cmed0,
                     inp.dhe,inp.cpsl,inp.rhosl,inp.spb,inp.spc,inp.ts,inp.qs,
                     inp.as,inp.tsd,inp.qtis,inp.hs,inp.tav,inp.xffm,inp.zp,
                     inp.z0,inp.za,inp.ua,inp.ta,inp.rh,inp.stab,inp.ala;
                     solver=solver, steady_solver_kwargs=steady_solver_kwargs,
                     kwargs...)
end

function slab_main(idspl::I,ncalc::I,wms::F,cps::F,tbp::F,cmed0::F,dhe::F,cpsl::F,rhosl::F,
                   spb::F,spc::F,ts::F,qs::F,as::F,tsd::F,qtis::F,hs::F,tav::F,xffm::F,
                   zp::AbstractVector{F},z0::F,za::F,ua::F,ta::F,rh::F,stab::F,
                   ala::F;msfm::I=11,mnfm::I=50,mffm::I=61,
                   solver=SLABLegacySolver(), steady_solver_kwargs=(;)) where {
                   I <: Integer, F <: AbstractFloat}

    #c  number of zp values
    # nzpm = 1
    # for i in 2:4
    #     if zp[i] == zero(F)
    #         break
    #     else
    #         nzpm = i
    #     end
    # end

   # select appropriate release type
    if idspl == 3
        # vertical jet
        idpf,nxtr,vecs,vars,params,dt = _slab_init_vjet(3,ncalc,msfm,mnfm,mffm,wms,cps,tbp,cmed0,
                                                    dhe,cpsl,rhosl,spb,spc,ts,qs,as,tsd,qtis,hs,tav,
                                                    xffm,zp,z0,za,ua,ta,rh,stab,ala)
        if idpf < 2
            phases = _slab_int_steady_state_impl!(vecs, vars, params, idpf, nxtr;
                                                  solver=solver, solver_kwargs=steady_solver_kwargs)
            steady_vecs = phases.steady_state
            steady_ode_solution = phases.ode_solution
            steady_ode_segments = phases.ode_segments
            if phases.transient_started
                _slab_int_transient!(vecs,phases.transient_vars,params,2,
                                     phases.nxtr,phases.transient_dt)
                transient_vecs = vecs
                transient_vars = phases.transient_vars
                transient_index = phases.nxtr
            else
                transient_vecs = nothing
                transient_vars = nothing
                transient_index = nxtr
            end
        else
            _slab_int_transient!(vecs,vars,params,idpf,nxtr,dt)
            steady_vecs = nothing
            steady_ode_solution = nothing
            steady_ode_segments = nothing
            transient_vecs = vecs
            transient_vars = vars
            transient_index = nxtr
        end
    else
        # default is a horizontal jet
        idpf,nxtr,vecs,vars,params,dt = _slab_init_hjet(2,ncalc,msfm,mnfm,mffm,wms,cps,tbp,cmed0,
                                                    dhe,cpsl,rhosl,spb,spc,ts,qs,as,tsd,qtis,hs,tav,
                                                    xffm,zp,z0,za,ua,ta,rh,stab,ala)
        if idpf < 2
            phases = _slab_int_steady_state_impl!(vecs, vars, params, idpf, nxtr;
                                                  solver=solver, solver_kwargs=steady_solver_kwargs)
            steady_vecs = phases.steady_state
            steady_ode_solution = phases.ode_solution
            steady_ode_segments = phases.ode_segments
            if phases.transient_started
                _slab_int_transient!(vecs,phases.transient_vars,params,2,
                                     phases.nxtr,phases.transient_dt)
                transient_vecs = vecs
                transient_vars = phases.transient_vars
                transient_index = phases.nxtr
            else
                transient_vecs = nothing
                transient_vars = nothing
                transient_index = nxtr
            end
        else
            _slab_int_transient!(vecs,vars,params,idpf,nxtr,dt)
            steady_vecs = nothing
            steady_ode_solution = nothing
            steady_ode_segments = nothing
            transient_vecs = vecs
            transient_vars = vars
            transient_index = nxtr
        end
    end

    cc_vecs = editcc(vecs,params,mffm)

    steady_solution = steady_vecs === nothing ? nothing :
        SLAB_Steady_Solution(params,steady_vecs,editcc(steady_vecs,params,mffm),
                             _slab_initial_steady_state(steady_vecs,vars),
                             steady_ode_solution,steady_ode_segments)
    transient_solution = transient_vecs === nothing ? nothing :
        SLAB_Transient_Solution(params,transient_vecs,cc_vecs,
            _slab_initial_transient_state(transient_vecs,transient_vars,
                clamp(transient_index,1,length(transient_vecs.x))))

    interpolations = _slab_output_interpolations(cc_vecs,steady_solution,params)

    return SLAB_Output(params,vecs,cc_vecs,interpolations,steady_solution,transient_solution)

end

end