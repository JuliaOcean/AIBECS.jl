using Documenter

# Ensure that the DataDeps download work remotely
ENV["JULIA_DEBUG"] = "Documenter"
ENV["DATADEPS_ALWAYS_ACCEPT"] = true
ENV["GKSwstype"] = "100"

deployconfig = Documenter.auto_detect_deploy_system()
Documenter.post_status(deployconfig; type="pending", repo="github.com/JuliaOcean/AIBECS.jl.git")

using Literate
using AIBECS

# generate tutorials and how-to guides using Literate
src = joinpath(@__DIR__, "src")
lit = joinpath(@__DIR__, "lit")
notebooks = joinpath(src, "notebooks")

execute = false # Set to true for executing notebooks and documenter!
nb = true      # Set to true to generate the notebooks
for (root, _, files) in walkdir(lit), file in files
    splitext(file)[2] == ".jl" || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, lit=>src))[1]
    Literate.markdown(ipath, opath, flavor = Literate.DocumenterFlavor())
    # Literate.script(ipath, opath, execute = execute)
    # nb && Literate.notebook(ipath, notebooks, execute = execute)
end

# Documentation structure
ismd(f) = splitext(f)[2] == ".md"
pages(folder) = [joinpath(folder, f) for f in readdir(joinpath(src, folder)) if ismd(f)]

makedocs(
    sitename="AIBECS.jl",
    doctest = false, # TODO guessing I should remove that when actually deploying?
    # options
    modules = [AIBECS],
    # organisation
    pages = Any[
        "Home" => "index.md",
        "Tutorials" => pages("tutorials"),
        "How-to guides" => pages("howtos"),
        "Explanation" => pages("explanation"),
        "Reference" => pages("reference")
        ],
    warnonly = true,
    format = Documenter.HTML(prettyurls = true),
    remotes = nothing
)

# Deploy
deploydocs(
    repo = "github.com/JuliaOcean/AIBECS.jl.git",
    push_preview = true
)


#=
To edit locally, make sure execute is set to false and run

using LiveServer
servedocs(literate=joinpath("docs", "lit"), doc_env=true)

from the root of the package in development
=#
