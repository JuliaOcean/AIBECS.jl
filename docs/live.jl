using LiveServer

sw = LiveServer.SimpleWatcher()

# source files to check for updates to avoid infinite loop
EXAMPLEDIR = joinpath(@__DIR__, "src", "examples")
EXAMPLES_jl = [f for f in readdir(EXAMPLEDIR) if endswith(f, ".jl")]
SOURCES = [
    abspath(joinpath(@__DIR__, "make.jl"))
    abspath(joinpath(@__DIR__, "src", "index.md"))
    abspath(joinpath(@__DIR__, "src", "prerequisites.md"))
    [abspath(joinpath(EXAMPLEDIR, example)) for example in EXAMPLES_jl]
]

append!(sw.watchedfiles, LiveServer.WatchedFile.(SOURCES))

function callback(x)
    # only trigger for source files to avoid infinite loop
    if x in SOURCES
        include(abspath(joinpath(@__DIR__, "make.jl")))
        LiveServer.file_changed_callback(x)
    end
end

callback(abspath(joinpath(@__DIR__, "make.jl"))) # make sure files exist
LiveServer.set_callback!(sw, callback)

LiveServer.serve(sw; dir = abspath(joinpath(@__DIR__, "build")), verbose=true)



