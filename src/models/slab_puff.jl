include("slab/slab.jl")

using .slab

# defining type for dispatch
struct SLAB <: PuffModel end

struct SLABSolution <: Puff
    out::SLAB_Output
end

function puff(scenario::Scenario, ::Type{SLAB}, eqs::EquationSet=DefaultSet();)

end