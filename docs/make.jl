using Documenter
using DocumenterVitepress

# Ensure that the DataDeps download work remotely
ENV["JULIA_DEBUG"] = "Documenter"
ENV["DATADEPS_ALWAYS_ACCEPT"] = true
ENV["GKSwstype"] = "100"

using Literate
using AIBECS
using Plots                 # activate AIBECSRecipesBaseExt so plot-recipe stubs gain methods
using Distributions  # activate parameter extensions for @docs introspection
using Bijectors      # activate parameter extensions for @docs introspection
using DataFrames     # activate parameter extensions for @docs introspection

# generate tutorials and how-to guides using Literate
src = joinpath(@__DIR__, "src")
lit = joinpath(@__DIR__, "lit")
notebooks = joinpath(src, "notebooks")

execute = false # Set to true for executing notebooks and documenter!
nb = true      # Set to true to generate the notebooks
for (root, _, files) in walkdir(lit), file in files
    splitext(file)[2] == ".jl" || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, lit => src))[1]
    Literate.markdown(ipath, opath, flavor = Literate.DocumenterFlavor())
    # Literate.script(ipath, opath, execute = execute)
    # nb && Literate.notebook(ipath, notebooks, execute = execute)
end

# Documentation structure
ismd(f) = splitext(f)[2] == ".md"
pages(folder) = [joinpath(folder, f) for f in readdir(joinpath(src, folder)) if ismd(f)]

makedocs(
    sitename = "AIBECS.jl",
    doctest = false, # TODO guessing I should remove that when actually deploying?
    draft = get(ENV, "DRAFT", "false") == "true", # DRAFT=true skips @example/@repl/@setup/@eval blocks for fast markdown iteration
    # options
    modules = [AIBECS],
    # organisation
    pages = Any[
        "Home" => "index.md",
        "Tutorials" => pages("tutorials"),
        "How-to guides" => pages("howtos"),
        "Explanation" => pages("explanation"),
        "Reference" => pages("reference"),
        "Publications" => pages("publications"),
    ],
    warnonly = [:missing_docs],   # internals are intentionally omitted from the curated reference page
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "https://github.com/JuliaOcean/AIBECS.jl",
        devbranch = "main",
        devurl = "dev",
        build_vitepress = get(ENV, "CI", "false") == "true",
    ),
    remotes = nothing
)

# Deploy (only on CI — locally there is nothing to deploy to)
if get(ENV, "CI", "false") == "true"
    DocumenterVitepress.deploydocs(
        repo = "github.com/JuliaOcean/AIBECS.jl.git",
        target = joinpath(@__DIR__, "build"),
        branch = "gh-pages",
        devbranch = "main",
        push_preview = true,
    )
end

#=
To render and preview the docs locally:

Recommended — VitePress dev server with hot reload (DocumenterVitepress-native):
    1. julia --project=docs -e 'include("docs/make.jl")'    # generates VitePress sources
       (locally `build_vitepress` is automatically false — see the format above —
       so this step skips the production build and just emits the sources.)
    2. julia --project=docs -e 'using DocumenterVitepress; DocumenterVitepress.dev_docs("docs/build/1")'
    Browser opens at http://localhost:5173 with hot module reload.
    On any src/ edit: re-run step 1; the browser auto-refreshes.

Note: since DocumenterVitepress 0.2, the local build lands in
`docs/build/1/` (numbered subfolders, one per base — for local dev there
is only one, with base = ""). That's why `dev_docs` points at `build/1`
rather than `build/`.

Fast markdown-only iteration (skip @example/@repl/@setup/@eval evaluation):
    DRAFT=true julia --project=docs -e 'include("docs/make.jl")'
    Code blocks render as plain code — useful when you only want to check
    page structure, prose, layout, or link targets. Drop `DRAFT=true` to
    evaluate code blocks again.

Note on `LiveServer.servedocs`: it auto-rebuilds on src/ changes but cannot
preview a VitePress site correctly (asset paths are absolute under
/AIBECS.jl/, and servedocs serves docs/build/ rather than the rendered
docs/build/final_site/). Use `dev_docs` instead for visual preview.
=#
