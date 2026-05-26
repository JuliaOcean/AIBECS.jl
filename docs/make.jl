using Documenter
using DocumenterVitepress

# Ensure that the DataDeps download work remotely
ENV["JULIA_DEBUG"] = "Documenter"
ENV["DATADEPS_ALWAYS_ACCEPT"] = true
ENV["GKSwstype"] = "100"

using Literate
using AIBECS
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

# Pre-compute the inline benchmark-table markdown so `explanation/solvers.md`'s
# `@eval` block can embed it. Documenter changes cwd / `@__DIR__` to the build
# directory when evaluating `@eval` blocks, so reading the source file from
# inside the block fails on a fresh build. Read it here, then:
#   1. strip the autogen HTML comment,
#   2. drop the "# Tier: small" metadata header (the surrounding section in
#      solvers.md already covers it),
#   3. demote the per-circulation `##` headings to bold lead-in paragraphs.
# Step 3 is a workaround for DocumenterVitepress: it flattens parsed-MD
# heading levels to H1 when emitting the rendered .md, which would
# collide with the page-title H1. Bold paragraphs keep the visual
# hierarchy without the level collision.
const __LATEST_TIMINGS_SMALL_MD = let
    raw = read(joinpath(src, "explanation", "latest_timings_small.md"), String)
    raw = replace(raw, r"<!--[^\n]*-->\n?" => "")
    raw = replace(raw, r"^# Tier:[^\n]*\n+"m => "")
    raw = replace(raw, r"^## (.+)$"m => s"**\1**")
    raw
end

makedocs(
    sitename = "AIBECS.jl",
    doctest = false, # TODO guessing I should remove that when actually deploying?
    draft = get(ENV, "DRAFT", "false") == "true", # DRAFT=true skips @example/@repl/@setup/@eval blocks for fast markdown iteration
    # options
    modules = [AIBECS],
    # organisation
    pages = Any[
        "Home" => "index.md",
        "Tutorials" => [
            "tutorials/1_ideal_age.md",
            "tutorials/2_radiocarbon.md",
            "tutorials/3_Pmodel.md",
            "tutorials/4_dustmodel.md",
            "tutorials/5_river_discharge.md",
            "tutorials/6_groundwater_discharge.md",
        ],
        "How-to guides" => [
            "howtos/1_parameters.md",
            "howtos/2_plot.md",
            "howtos/2b_plot_makie.md",
            "howtos/3_cruiseplot.md",
            "howtos/4_parameter_optimization.md",
            "howtos/5_fluxes.md",
            "howtos/6_sinking_particles.md",
            "howtos/7_etopo.md",
            "howtos/8_nonlinearsolve.md",
            "howtos/analytical_derivatives.md",
        ],
        "Explanation" => [
            "explanation/1_concept.md",
            "explanation/2_tracer_transport_operators.md",
            "explanation/solvers.md",
            "explanation/datasets.md",
        ],
        "Reference" => ["reference/functions.md"],
        "Roadmap" => ["roadmap.md"],
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
