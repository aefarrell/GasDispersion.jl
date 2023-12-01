include("slab/slab.jl")

using .slab

# defining type for dispatch
struct SLAB <: PuffModel end

struct SLABSolution{I <: Integer, F <: Number, A <: AbstractVector{F}} <: Puff
    out::SLAB_Output{I,F,A}
end

function puff(scenario::Scenario, ::Type{SLAB}, eqs::EquationSet=DefaultSet();)

end