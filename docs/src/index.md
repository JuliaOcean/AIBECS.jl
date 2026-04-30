```@raw html
---
layout: home

hero:
  name: AIBECS.jl
  text:
  tagline: Algebraic Implicit Biogeochemistry Elemental Cycling System
  image:
    src: /assets/logo.png
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

!!! note
    The AIBECS is being developed primarily by Benoît Pasquier (postdoc with Seth John at the University of Southern California).
    If you use the AIBECS in your research, [please cite it](https://doi.org/10.5281/zenodo.2864051)!
    Similarly, if you access data from within AIBECS (like the OCIM or OCCA ocean circulations) please cite them too.

!!! warning
    This package is in active development, so nothing is set in stone, and things may be broken sometimes.
    And if you have any suggestions or feature requests, do not hesitate to [start an issue](https://github.com/JuliaOcean/AIBECS.jl/issues), or better even, [submit a pull request](https://github.com/JuliaOcean/AIBECS.jl/pulls)!
