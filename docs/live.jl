using LiveServer

sw = LiveServer.SimpleWatcher()

# source files to check for updates to avoid infinite loop
EXAMPLEDIR = joinpath(@__DIR__, "src", "examples")
EXAMPLES_jl = [f for f in readdir(EXAMPLEDIR) if endswith(f, ".jl")]
SOURCES = [
    joinpath(@__DIR__, "make.jl")
    joinpath(@__DIR__, "src", "index.md")
    joinpath(@__DIR__, "src", "prerequisites.md")
    [joinpath(EXAMPLEDIR, example) for example in EXAMPLES_jl]
]

append!(sw.watchedfiles, LiveServer.WatchedFile.(SOURCES))

function callback(x)
    # only trigger for source files to avoid infinite loop
    if x in SOURCES
        include(joinpath(@__DIR__, "make.jl"))
    end
end

callback(joinpath(@__DIR__, "make.jl")) # make sure files exist
LiveServer.set_callback!(sw, callback)

LiveServer.serve(sw; dir = joinpath(@__DIR__, "build"), verbose=true)



