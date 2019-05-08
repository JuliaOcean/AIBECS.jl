using Documenter, AIBECS

makedocs(
    sitename="AIBECS Documentation",
    # options
    modules = [AIBECS]
)

deploydocs(
    repo = "github.com/briochemc/AIBECS.jl.git",
)