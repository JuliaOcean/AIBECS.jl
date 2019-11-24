using Documenter, AIBECS
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

# Generate examples
include("generate.jl")

tutorials_jl = [f for f in readdir(tutorials_directory) if endswith(f, ".jl")]
howtos_jl = [f for f in readdir(howtos_directory) if endswith(f, ".jl")]

tutorials_md = [replace(f, ".jl" => ".md") for f in tutorials_jl]
howtos_md = [replace(f, ".jl" => ".md") for f in howtos_jl]

generated_tutorials = [joinpath("tutorials", "generated", f) for f in tutorials_md]
generated_howtos = [joinpath("howtos", "generated", f) for f in howtos_md]

makedocs(
    sitename="AIBECS.jl",
    doctest = false, # TODO guessing I should remove that when actually deploying?
    # options
    modules = [AIBECS],
    # organisation
    pages = Any[
        "Home" => "index.md",
        "Tutorials" => generated_tutorials,
        "How-to guides" => generated_howtos,
        "Explanation" => "explanation.md",
        "Reference" => "reference.md"
        ]
)

# Deploy
deploydocs(
    repo = "github.com/briochemc/AIBECS.jl.git",
)

