using Documenter, AIBECS

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
        "Examples" => GENERATEDEXAMPLES
        ]
)

# make sure there are no *.vtu files left around from the build
# Still do not know what these are
cd(joinpath(@__DIR__, "build", "examples", "generated")) do
    foreach(file -> endswith(file, ".vtu") && rm(file), readdir())
end

# Deploy
deploydocs(
    repo = "github.com/briochemc/AIBECS.jl.git",
)

