# generate tutorials and how-to guides
using Literate

# Directory for source files (in Literate format) of tutorials and how-to guides
tutorials_directory = abspath(joinpath(@__DIR__, "src", "tutorials"))
howtos_directory = abspath(joinpath(@__DIR__, "src", "howtos"))

# function to generate the files from the source directory
function generate(source_dir; execute = false)
    for source in readdir(source_dir)
        ~endswith(source, ".jl") && continue
        input = joinpath(source_dir, source)
        generated_dir = joinpath(source_dir, "generated")
        script = Literate.script(input, generated_dir)
        code = strip(read(script, String))
        if execute
            mdpost(str) = replace(str, "@__CODE__" => code)
            Literate.notebook(input, generated_dir, execute = true)
            Literate.markdown(input, generated_dir, postprocess = mdpost)
        else
            Literate.notebook(input, generated_dir, execute = false)
            Literate.markdown(input, generated_dir, documenter = false)
        end
    end
end
# Generate tutorials
generate(tutorials_directory; execute=false)

# Generate how-to guides
#generate(howtos_directory; execute=false)

