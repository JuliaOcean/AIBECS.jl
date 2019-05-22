## Prerequisites

This section provides a set of instructions and prerequisites to use AIBECS.

### 1. Install Julia

First things first, you must install [Julia](https://julialang.org). Click on the [Julia](https://julialang.org) link, look for the "download" buttons, and install the correct version for your OS.
Once this is done, you should be able to start Julia by typing

```bash
julia
```

in the terminal.
If not, find the Julia executable, and simply double click on it!
This should open a terminal session, and display something like this:

```julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.1.1 (2019-05-16)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia>
```

This is called the Julia REPL (for Read Eval Loop Print) and is used for interactive use of Julia.
Great job, Julia is now running on your computer! Congratulations!

If you want to learn more about Julia, you can read [the documentation](https://docs.julialang.org/en/v1/), there is a [Discourse forum](https://discourse.julialang.org/), and there is a [Slack channel](https://julialang.slack.com/messages) if you need help.
But for now you should not need any of those: The notebook will just require you to press Shift + Enter a couple of times.

### 2. Julia packages required for the notebooks

In Julia, you can access the package manager by simply typing `]` in the REPL.
Once you type `]`, the REPL changes to

```julia
(v1.1) pkg>
```

This means you're in the package-manager (or `pkg`) mode.

Note that you can exit the `pkg` mode by pressing the `delete` key, and this will revert the Julia prompt to its original form:

```julia
julia>
```

The packages you should install are:

- **AIBECS** (mandatory)

    To create a global steady-state biogeochemistry model, we will be using the [AIBECS](https://github.com/briochemc/AIBECS.jl) package.
    You install it, via

    ```julia
    add AIBECS
    ```

    in `pkg` mode, which should look like

    ```julia
    (v1.1) pkg> add AIBECS
    ```

    This should only take a few seconds.

- If you want nice-looking maps, install **Cartopy**

    In order to plot things, i.e., to look at the output of the beautiful work you will be doing with AIBECS, you will need a plotting package.
    I suggest using Python's [Cartopy](https://scitools.org.uk/cartopy/docs/latest/) because, well, it looks pretty.
    You can install it via the `Conda` package from within Julia.
    Install it via `add Conda` in `pkg` mode.
    You should see something like:

    ```julia
    (v1.1) pkg> add Conda
    ```

    Then back to the normal REPL mode (press `escape` to get the normal Julia REPL prompt), install cartopy within conda (within Julia) via

    ```julia
    julia> using Conda

    julia> Conda.add("Cartopy")
    ```


    Then Install PyPlot:
    ```julia
    (v1.1) pkg> add PyPlot
    ```

    And PyCall:
    ```julia
    (v1.1) pkg> add PyCall
    ```

    This should only take a few seconds as well.

    !!! note
        You may want to install Cartopy differently, or even use a different plotting package.
        This is merely a suggestion that has worked well for me.

- If you want to run the notebooks, install **IJulia**

    The [IJulia](https://github.com/JuliaLang/IJulia.jl) package allows you to launch JupyterLab from Julia.
    To install it, in `pkg` mode, type `add IJulia` (and press return), and you should see something like:

    ```julia
    (v1.1) pkg> add IJulia
    ```

    This should only take a few seconds as well.

    !!! note
        You could install JupyterLab externally to Julia and use that instead.
        The solution proposed here is to facilitate usage for those who do not already have or know how to install JupyterLab.

    !!! note
        You may need to build **CodecZlib** to run the notebooks.
        If you see an error mentioning you should build it, then build it!
        It's easy, just go in `pkg` mode and type `build CodecZlib`.
        It should look like this:
        ```julia
        (v1.1) pkg> build CodecZlib
        ```

If you followed all these steps you should be able to use the notebooks!

The final step is to start JupyterLab.
First, make sure you are in the normal Julia REPL mode (i.e., press `delete` if you are in `pkg` mode.)
Then, tell Julia that you want to "use" IJulia:

```julia
julia> using IJulia
```

!!! note
    You can just copy paste the code above (including the `julia>` bits), and the REPL will know to not paste those automatically.
    Everytime a package is used for the first time, Julia will precompile it (which can take a few seconds to minutes, depending on the package — don't worry, just let it finish).

Finally, you can start JupyterLab from Julia by simply typing `jupyerlab()` in Julia.
It should look like:

```julia
julia> jupyterlab()
```

!!! note
    If Julia asks you if you want Conda to install JupyterLab, just say "yes" (i.e., type `y`).
    After a couple seconds/minutes of downloads and installations, you should be all set up and a browser window should open with JupyterLab!

Just navigate to the notebook of your choice with JupyterLab in your browser and double-click on the notebook!


