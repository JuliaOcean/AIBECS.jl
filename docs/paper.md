---
title: 'AIBECS.jl: The ideal tool for exploring global marine biogeochemical cycles.'
tags:
  - Julia
  - Oceanography
  - Biogeochemical Cycles
  - Earth Sciences
authors:
  - name: Benoit Pasquier
    orcid: 0000-0002-3838-5976
    affiliation: "1" # (Multiple affiliations must be quoted)
affiliations:
 - name: Department of Earth Sciences, University of Southern California
   index: 1
date: 27 July 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

<!--
Has a clear description of the high-level functionality and purpose of the software for a diverse, non-specialist audience been provided?
--->

Our Ocean is a beautifully complex system that continuously transports and mixes elements on a multitude of timescales and lengthscales.
In order to improve our undesrstanding of global biogeochemical cycles, oceanographers collect more data every day that are then used in an increasingly wide variety of numerical models.
While the simplest "toy" box models that can be computed by hand are fairly accessible, the most advanced fine-resolution three-dimensonal models often require high-performance computing clusters and computational-science expertise.
A high barrier to entry for modelling in oceanography aggravates the workload of novices, digs a wedge between sea-going oceanographers and modellers, and hinders the advances in the field.
[AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl) [@Pasquier:2021] is a [JuliaOcean](https://github.com/JuliaOcean/)-affiliated package written in [Julia](https://julialang.org/) [@Bezanson_etal:2017], which aims to lower this barrier by providing an easy-to-use, open-source, and composable framework for simulating global marine tracers.






A conceptual model of the cycle of any marine tracer essentially requires two components.
(i) A model of how the tracer is transported (be it the ocean currents and eddies, gravitational settling, or a combination of those), and (ii) a model of the local external inputs and outputs at any location.
AIBECS.jl is built on this concept and allows users to build numerical models of marine tracers by selecting a circulation (and/or vertical transport) and local sources and sinks.
Tools for generating the equations, solving them, and then diagnosing and plotting the simulated tracers are also provided, either directly by AIBECS, its dependencies, or by satellite packages in the AIBECS and Julia ecosystem.


To the best of my knowledge, only the excellent [AWESOME OCIM](https://github.com/MTEL-USC/AWESOME-OCIM) [A Working Environment for Simulating Ocean Movement and Elemental cycling within the Ocean Circulation Inverse Model, @John_etal_2020] provides a comparable framework.
However, the AWESOME OCIM is written in a proprietary language (MATLAB), making the AIBECS the only truly open-source package for simulating marine tracers easily and for free.
AIBECS also provides a number of advantages over its predecessor, including the ability to simulate mutiple tracers and nonlinear systems, integrated unit-conversions, and plug-ins for autodifferentiation and optimization, to name a few.
Simple global biogeochemical models can be built in a few lines of code, making AIBECS.jl a useful tool for teaching and for exploratory reserach.
Its more advanced capabilities also allows more experienced coders to use AIBECS for cutting-edge research.



Through a simple user interface, AIBECS.jl provides access to a variety of steady-state ocean circulation models.
These currently include the Ocean Circulation Inverse Model (OCIM) v0.1 [@Primeau_etal:2013], v1.0 [@DeVries_Primeau:2011; @DeVries:2014], and v2.0 [@DeVries_Holzer:2019], the MITgcm-built Ocean Comprehensible Atlas (OCCA) ocean-state estimate model [@Forget:2010], of which the download is handled by the [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl) package [@DataDeps:2018; @White_etal:2019].
AIBECS.jl also offers classic two-box and three-box models [@Sarmiento_Gruber:2006; @Archer_etal:2000].
The [OceanGrids.jl](https://github.com/briochemc/OceanGrids.jl) package, on which AIBECS depends, provides the underlying grid configuration types as well as regridding and interpolating routines.
Licensing and authorizations permitting, additional circulation models will be added to the collection.
Swapping the underlying circulation model and grid requires a single-line-of-code change, facilitating benchmarks among them and, given past- or future-ocean circulation models, could be used in paleoceanography and climate-change studies.



AIBECS.jl also provides extra functionality to facilitate the generation of numerical models.
Tooling to simulate gravitational settling of non-buoyant particulate tracers is provided by the `transportoperator` function.
In complement, AIBECS provides access to a number of predefined fields that can be used to generate source and sink processes.
Fine-resolution (1-arc-minute) topography from the ETOPO1 dataset [@Amante_Eakins_2009] can be used for a refined interception of particulate fluxes by subgrid topographic features not captured by coarser circulation models.
For aeolian deposition, AIBECS.jl includes aerosol-type- and region-of-origin-partitioned dust deposition fields [@Chien_etal_2016; @Kok_etal_2021b].
Datasets for global river discharge [@Dai_2017; @Dai_Trenberth_2002] and surface groundwater discharge [@Luijendijk_etal_2019; @Luijendijk_etal_2020] are included.
For hydrothermal-sourced tracers, the helium fluxes from the Earth's mantle computed with the OCIM v1.0 and v2.0 are available when loading the corresponding circulation models [@DeVries:2014; DeVries_Holzer:2019].
AIBECS.jl also provides access to the data included with the AWESOME OCIM framework [@John_etal_2020], namely data from the Global Ocean Data Analysis Project [GLODAP, @Lauvset_etal_2016; @Olsen_etal_2016], P-cycling modelled fields from [@Weber_etal_Science_2018], nepholoid layers [@Gardner_etal_2018a; @Gardner_etal_2018b; Taburet_etal_2019], as well as other data already present within AIBECS or sattelite packages.
Also useful to global biogeochemistry modelling are data from the World Ocean Atlas [@WOA_2018_nut] that can be downloaded, assisted by external package [WorldOceanAtlasTools.jl](https://github.com/briochemc/WorldOceanAtlasTools.jl) [@WorldOceanAtlasTools.jl-2019].
Similarly GEOTRACES data [@Schlitzer:2018df] can be handled by the [GEOTRACES.jl](https://github.com/briochemc/GEOTRACES.jl) package (although GEOTRACES requires manual download of the data).
Finally, more advanced usage such as optimization is facilitated by the [F1Method.jl](https://github.com/briochemc/F1Method.jl) package [@F1Method], which provides efficient gradient and Hessian computations of objective functions defined through AIBECS.jl, which can then be directly fed to optimization routines such as, e.g., the [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) package [@Optim.jl-2018] to optimize biogeochemical parameters.






Internally, AIBECS.jl uses a quasi-Newton solver [@Kelley_2003_1] translated from MATLAB to Julia and taylored to the context of marine tracers to solve/simulate tracers.
AIBECS uses forward-mode auto-differentiation from the [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) package [@RevelsLubinPapamarkou2016] for the nonlinear parts of the PDE to generate the Jacobian required for the solver.
In future versions, AIBECS will target wrapper packages, such as [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl) from the [SciML](https://github.com/SciML) organization, to be able to plug into other nonlinear solvers available in the Julia ecosystem.





Model parameters are currently handled by a combination of dependencies, which include [UnPack.jl](https://github.com/mauro3/Unpack.jl), [FieldMetadata.jl](https://github.com/rafaqz/FieldMetadata.jl), [Flatten.jl](https://github.com/rafaqz/Flatten.jl), [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl), [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) [@Besancon_etal_2021; @Distributions.jl-2019], [Unitful.jl](https://github.com/PainterQubits/Unitful.jl), and [Bijectors.jl](https://github.com/TuringLang/Bijectors.jl).
These dependencies allow for parameters with attached metadata such as units and prior distributions to be used within the AIBECS.jl framework.






The AIBECS.jl package is registered with Julia's [default package registry (GENERAL)](https://github.com/JuliaRegistries/General), such that installation takes a single line of code from within Julia.
The package documentation, which is built through continuous integration (CI), includes runnable tutorial examples that are both available online for consultation and as Jupyter notebooks.
Continuous integration through GitHub actions also includes a fairly complete suite of tests.







# Statement of need


<!--
Do the authors clearly state what problems the software is designed to solve and who the target audience is?

Do the authors describe how this software compares to other commonly-used packages?

Is the paper well written (i.e., it does not require editing for structure, language, or writing quality)?

References: Is the list of references complete, and is everything cited appropriately that should be cited (e.g., papers, datasets, software)? Do references in the text use the proper citation syntax?
-->




Steady-state ocean circulation models, which can bypass the need for long and costly spinups altogether, have been increasingly used because they provide significant computational efficiency gains [@Kwon_Primeau_2006; @Kwon_Primeau_2008; @DeVries_etal_GBC_2013; @Frants_etal_2016; @DeVries_Weber_GBC_2017; @Pasquier_Holzer_2017].

<!--
Need to find how citaations are ordered for prefix to work
-->

 [Primeau, 2003; Khatiwala et al., 2005; Primeau, 2005, Kwon and Primeau, 2006; Khatiwala et al., 2009; Weber and Deutsch, 2010; Wang et al., 2019], with existing and ongoing efforts to extract steady-state circulation estimates from existing general circulation models [Bardin, Kvale et al., 2017, Chamberlain et al., 2019, Zanna et al., 2019, competition for unswjob]

AO stuff:
CO2 uptake (DeVries, 2014; DeVries et al., 2017), reveal spatial variability in phytoplankton carbon to nutrient uptake ratios (DeVries and Deutsch, 2014; Teng et al., 2014), and constrain the export and fate of organic carbon (Weber et al., 2016; Roshan and DeVries, 2017; DeVries and Weber, 2017).
The OCIM has also been used to constrain the global cycling and distribution of a wide range of elements including phosphorous (Primeau et al., 2013; DeVries et al., 2014), silicon (Holzer et al., 2014; DeVries et al., 2017), nitrogen (DeVries et al., 2012; Weber and Deutsch, 2012; DeVries et al., 2013; Weber and Deutsch, 2014; Wang
et al., 2019), and zinc (Roshan et al., 2018; Weber et al., 2018)

AIBECS.jl aims to facilitate and standardize the use of these models by providing
(i) integrated management of ocean-circulation-model data,
(ii) handling regridding and inteprolation from model one ocean-model grid to another,
(iii) ...

OCCA
OCIM
AWESOME OCIM



<!---
`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.
--->

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References