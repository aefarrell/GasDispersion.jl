# original SLAB subroutines

include("subroutines/slope.jl")
include("subroutines/slopepf.jl")
include("subroutines/solve.jl")
include("subroutines/solvepf.jl")
include("subroutines/thermo.jl")
include("subroutines/eval.jl")
include("subroutines/evalpf.jl")
include("subroutines/entran.jl")
include("subroutines/store.jl")
include("subroutines/editcc.jl")

# integration routines
include("subroutines/integrate_steadystate.jl")
include("subroutines/integrate_transient.jl")

# intialize parameters
include("subroutines/initialize_common.jl")
include("subroutines/initialize_hjet.jl")
include("subroutines/initialize_vjet.jl")