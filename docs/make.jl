using Documenter, Literate, AIBECS

# Ensure that the DataDeps download work remotely
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

# generate tutorials and how-to guides using Literate
src = joinpath(@__DIR__, "src")
lit = joinpath(@__DIR__, "lit")
notebooks = joinpath(@__DIR__, "notebooks")

execute = true # Set to true for executing notebooks and documenter!
nb = true      # Set to true to generate the notebooks
for (root, _, files) in walkdir(lit), file in files
    splitext(file)[2] == ".jl" || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, lit=>src))[1]
    Literate.markdown(ipath, opath, documenter = execute)
    if nb
        opath_nb = splitdir(replace(ipath, lit=>notebooks))[1]
        Literate.notebook(ipath, opath_nb, execute = execute)
    end
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
        ]
)

# Deploy
deploydocs(
    repo = "github.com/briochemc/AIBECS.jl.git",
)


#=
To edit locally, make sure execute is set to false and run 


using LiveServer
servedocs(literate=joinpath("docs", "lit"))

from the root of the package in development
=#
