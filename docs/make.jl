using Documenter, AIBECS
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

# Generate examples
include("generate.jl")

EXAMPLES_jl = [f for f in readdir(EXAMPLEDIR) if endswith(f, ".jl")]
EXAMPLES_md = [replace(f, ".jl" => ".md") for f in EXAMPLES_jl]
GENERATEDEXAMPLES = [joinpath("examples", "generated", f) for f in EXAMPLES_md]

makedocs(
    sitename="AIBECS.jl",
    doctest = false, # TODO guessing I should remove that when actually deploying?
    # options
    modules = [AIBECS],
    # organisation
    pages = Any[
        "Home" => "index.md",
        "Prerequisites" => "prerequisites.md",
        "Examples" => GENERATEDEXAMPLES,
        "Function index" => "functions.md"
        ]
)

# Deploy
deploydocs(
    repo = "github.com/briochemc/AIBECS.jl.git",
)

