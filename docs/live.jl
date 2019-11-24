using LiveServer

sw = LiveServer.SimpleWatcher()

# source files to check for updates to avoid infinite loop
tutorials_directory = joinpath(@__DIR__, "src", "tutorials")
tutorials_jl = [f for f in readdir(tutorials_directory) if endswith(f, ".jl")]
howtos_directory = joinpath(@__DIR__, "src", "howtos")
howtos_jl = [f for f in readdir(howtos_directory) if endswith(f, ".jl")]
sources = [
    abspath(joinpath(@__DIR__, "make.jl"))
    abspath(joinpath(@__DIR__, "src", "index.md"))
    abspath(joinpath(@__DIR__, "src", "prerequisites.md"))
    abspath(joinpath(@__DIR__, "src", "functions.md"))
    [abspath(joinpath(tutorials_directory, example)) for example in tutorials_jl]
    [abspath(joinpath(howtos_directory, example)) for example in howtos_jl]
]

append!(sw.watchedfiles, LiveServer.WatchedFile.(sources))

function callback(x)
    # only trigger for source files to avoid infinite loop
    if x in sources
        include(abspath(joinpath(@__DIR__, "make.jl")))
        LiveServer.file_changed_callback(x)
    end
end

callback(abspath(joinpath(@__DIR__, "make.jl"))) # make sure files exist
LiveServer.set_callback!(sw, callback)

LiveServer.serve(sw; dir = abspath(joinpath(@__DIR__, "build")), verbose=true)



