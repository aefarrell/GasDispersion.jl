push!(LOAD_PATH,"../src/")

using Documenter, GasDispersion

ENV["GKSwstype"] = "100"

makedocs(sitename="GasDispersion.jl Documentation",
         pages=["scenarios.md", "plume.md", "puff.md", "equation_sets.md", "references.md", "function_index.md"])

deploydocs(
    repo = "github.com/aefarrell/GasDispersion.jl.git",
)
