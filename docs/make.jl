push!(LOAD_PATH,"../src/")

using Documenter, GasDispersion

makedocs(sitename="GasDispersion.jl Documentation")

deploydocs(
    repo = "github.com/aefarrell/GasDispersion.jl.git",
)
