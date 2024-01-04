push!(LOAD_PATH,"../src/")

using Documenter, GasDispersion

makedocs(sitename="GasDispersion.jl Documentation",
         pages=["scenarios.md", "plume.md", "puff.md", "equation_sets.md", "function_index.md"])

deploydocs(
    repo = "github.com/aefarrell/GasDispersion.jl.git",
)
