This ReadMe is **not** part of the AIBECS documentation.
It is merely for myself to remember how to edit the documentation.

To edit the documentation, it is better to avoid running the full suite of scripts and examples.
In the `generate.jl` script, comment/uncomment the Literate functions that generate the markdown and` notebooks files to not execute them.

Then, from the root of AIBECS, run
```bash
julia --project=docs --color=yes -L docs/live.jl
```
which should start a local website at `http://localhost:8000/` that will host the documentation and refresh on any save of the examples or other documentation source files.

And don't forget to turn the execution on for a last run before you push!
