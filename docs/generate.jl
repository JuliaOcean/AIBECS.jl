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
    #Literate.notebook(input, GENERATEDDIR, execute = true)
    #Literate.markdown(input, GENERATEDDIR, postprocess = mdpost)
    # Altenrate comments to avoid running the code when editing the text only!
    #Literate.notebook(input, GENERATEDDIR, execute = false)
    Literate.markdown(input, GENERATEDDIR, documenter = false)
end

