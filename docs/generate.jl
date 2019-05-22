# generate examples
using Literate

# Directory where the examples in Literate.jl format are
EXAMPLEDIR = joinpath(@__DIR__, "src", "examples")
# Directory where the examples will be generated
GENERATEDDIR = joinpath(@__DIR__, "src", "examples", "generated")

for example in readdir(EXAMPLEDIR)
    ~endswith(example, ".jl") && continue
    input = abspath(joinpath(EXAMPLEDIR, example))
    script = Literate.script(input, GENERATEDDIR)
    code = strip(read(script, String))
    mdpost(str) = replace(str, "@__CODE__" => code)
    Literate.markdown(input, GENERATEDDIR, postprocess = mdpost)
    Literate.notebook(input, GENERATEDDIR, execute = true)
    # Altenrate comments to avoid running the code when editing the text only!
    #Literate.notebook(input, GENERATEDDIR, execute = false)
    #Literate.markdown(input, GENERATEDDIR, documenter=false)
end

# copy some figures to the build directory (no figures for me right now)
# cp(joinpath(@__DIR__, "../examples/figures/heat_square.png"),
#    joinpath(@__DIR__, "src/examples/generated/heat_equation.png"); force = true)
# 
# cp(joinpath(@__DIR__, "../examples/figures/coloring.png"),
#    joinpath(@__DIR__, "src/examples/generated/coloring.png"); force = true)

# remove any .vtu files in the generated dir (should not be deployed)
# I don't even know what VTU files are... Leaving this here for now.
cd(GENERATEDDIR) do
    foreach(file -> endswith(file, ".vtu") && rm(file), readdir())
end
