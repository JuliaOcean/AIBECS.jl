```@raw html
---
layout: home

hero:
  name: AIBECS.jl
  text:
  tagline: Algebraic Implicit Biogeochemistry Elemental Cycling System
  image:
    src: /logo.png
    alt: AIBECS logo
  actions:
    - theme: brand
      text: Get Started
      link: /tutorials/1_ideal_age
    - theme: alt
      text: View on GitHub
      link: https://github.com/JuliaOcean/AIBECS.jl

features:
  - icon: 📖
    title: Tutorials
    details: Step-by-step introductions to AIBECS, from ideal age to coupled models.
    link: /tutorials/1_ideal_age
  - icon: 🧭
    title: How-to Guides
    details: Goal-oriented walk-throughs for specific tasks like plotting and fluxes.
    link: /howtos/1_parameters
  - icon: 💡
    title: Explanation
    details: In-depth discussion of the concepts behind AIBECS.
    link: /explanation/1_concept
  - icon: 📚
    title: Reference
    details: Function docstrings and API reference.
    link: /reference/functions
---
```

!!! note
    If you are using the AIBECS for the first time, you must add it to your Julia environment, by typing
    ```
    ]add AIBECS
    ```
    at the REPL.

----

!!! note "Maintenance status"
    AIBECS is essentially in maintenance-only mode, and might remain so for a little while.
    Although it works (the software has been published and has been used in publications), I consider AIBECS still in development, with a bunch of features and fixes that are currently missing, but I do not have enough time to work on it unfortunately.
    Please get in touch if you have some bandwidth and/or funding and are interested in collaborating or developing AIBECS!

    If you use AIBECS in your research, [please cite it](https://doi.org/10.5281/zenodo.2864051)!
    Similarly, if you access data from within AIBECS (like the OCIM or OCCA ocean circulations) please cite them too.

!!! warning
    This package is in active development, so nothing is set in stone, and things may be broken sometimes.
    And if you have any suggestions or feature requests, do not hesitate to [start an issue](https://github.com/JuliaOcean/AIBECS.jl/issues), or better even, [submit a pull request](https://github.com/JuliaOcean/AIBECS.jl/pulls)!
