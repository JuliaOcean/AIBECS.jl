---
title: 'AIBECS.jl: The ideal tool for exploring global marine biogeochemical cycles.'
tags:
  - Julia
  - Oceanography
  - Biogeochemical Cycles
  - Earth Sciences
authors:
  - name: Benoît Pasquier^[corresponding author]
    orcid: 0000-0002-3838-5976
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Seth John
    orcid: 0000-0002-8257-626X
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: François Primeau
    orcid: 0000-0001-7452-9415
    affiliation: "2"
affiliations:
  - name: Department of Earth Sciences, University of Southern California
    index: 1
  - name: Department of Earth System Science, University of California, Irvine
    index: 2
date: 19 Aug 2021
bibliography: paper.bib

---

# Summary



Our Ocean is a beautifully complex system that continuously transports and mixes elements on a multitude of timescales and lengthscales.
In order to improve our understanding of global biogeochemical cycles, oceanographers collect more data every day that are then used in an increasingly wide variety of numerical models.
While the simplest "toy" box models that can be computed by hand are fairly accessible, the most advanced fine-resolution three-dimensional models often require high-performance computing clusters and computational-science expertise.
A high barrier to entry for modelling in oceanography aggravates the workload of novices, digs a wedge between sea-going oceanographers and modellers, and hinders the advances in the field.
[AIBECS.jl](https://github.com/JuliaOcean/AIBECS.jl) [@Pasquier:2021] is a [JuliaOcean](https://github.com/JuliaOcean/)-affiliated package written in [Julia](https://julialang.org/) [@Bezanson_etal:2017], which aims to lower this barrier by providing an easy-to-use, open-source, and composable framework for simulating global marine tracers.






A conceptual model of the cycle of any marine tracer essentially requires two components.
(i) A model of how the tracer is transported (be it the ocean currents and eddies, gravitational settling, or a combination of those), and (ii) a model of the local external inputs and outputs at any location.
AIBECS.jl is built on this concept and allows users to build numerical models of marine tracers by selecting a circulation (and/or vertical transport) and local sources and sinks.
Tools for generating the equations, solving them, and then diagnosing and plotting the simulated tracers are also provided, either directly by AIBECS, its dependencies, or by satellite packages in the AIBECS and Julia ecosystem.


To the best of our knowledge, only the excellent [AWESOME OCIM](https://github.com/MTEL-USC/AWESOME-OCIM) [A Working Environment for Simulating Ocean Movement and Elemental cycling within the Ocean Circulation Inverse Model, @John_etal_2020] provides a comparable framework.
However, the AWESOME OCIM is written in a proprietary language (MATLAB), making the AIBECS the only truly open-source package for simulating marine tracers easily and for free.
AIBECS also provides a number of advantages over its predecessor, including the ability to simulate multiple tracers and nonlinear systems, integrated unit-conversions, and plug-ins for autodifferentiation and optimization, to name a few.
Simple global biogeochemical models can be built in a few lines of code, making AIBECS.jl a useful tool for teaching and for exploratory research.
Its more advanced capabilities also allow more experienced coders to use AIBECS for cutting-edge research.



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
For hydrothermal-sourced tracers, the helium fluxes from the Earth's mantle computed with the OCIM v1.0 and v2.0 are available when loading the corresponding circulation models [@DeVries:2014; @DeVries_Holzer:2019].
AIBECS.jl also provides access to the data included with the AWESOME OCIM framework [@John_etal_2020], namely data from the Global Ocean Data Analysis Project [GLODAP, @Lauvset_etal_2016; @Olsen_etal_2016], P-cycling modelled fields from @Weber_etal_Science_2018, nepheloid layers [@Gardner_etal_2018a; @Gardner_etal_2018b; Taburet_etal_2019], as well as other data already present within AIBECS or satellite packages.
Also useful to global biogeochemistry modelling are data from the World Ocean Atlas [@WOA_2018_nut] that can be downloaded, assisted by external package [WorldOceanAtlasTools.jl](https://github.com/briochemc/WorldOceanAtlasTools.jl) [@WorldOceanAtlasTools.jl-2019].
Similarly, GEOTRACES data [@Schlitzer:2018df] can be handled by the [GEOTRACES.jl](https://github.com/briochemc/GEOTRACES.jl) package (although GEOTRACES requires manual download of the data).
More advanced usage such as optimization is facilitated by the [F1Method.jl](https://github.com/briochemc/F1Method.jl) package [@F1Method], which provides efficient gradient and Hessian computations of objective functions defined through AIBECS.jl, which can then be directly fed to optimization routines such as, e.g., the [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) package [@Optim.jl-2018] to optimize biogeochemical parameters.
Finally, plotting recipes for the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) [@Plots.jl], are available.






Internally, AIBECS.jl uses a quasi-Newton solver [@Kelley_2003_1] translated from MATLAB to Julia and tailored to the context of marine tracers to solve/simulate tracers.
AIBECS uses forward-mode auto-differentiation from the [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) package [@RevelsLubinPapamarkou2016] for the nonlinear parts of the PDE to generate the Jacobian required for the solver.
In future versions, AIBECS will target wrapper packages, such as [LinearSolve.jl](https://github.com/SciML/LinearSolve.jl) from the [SciML](https://github.com/SciML) organization, to be able to plug into other nonlinear solvers available in the Julia ecosystem.





Model parameters are currently handled by a combination of dependencies, which include [UnPack.jl](https://github.com/mauro3/Unpack.jl), [FieldMetadata.jl](https://github.com/rafaqz/FieldMetadata.jl), [Flatten.jl](https://github.com/rafaqz/Flatten.jl), [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl), [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) [@Besancon_etal_2021; @Distributions.jl-2019], [Unitful.jl](https://github.com/PainterQubits/Unitful.jl), and [Bijectors.jl](https://github.com/TuringLang/Bijectors.jl).
These dependencies allow for parameters with attached metadata such as units and prior distributions to be used within the AIBECS.jl framework.






The AIBECS.jl package is registered with Julia's [default package registry (GENERAL)](https://github.com/JuliaRegistries/General), such that installation takes a single line of code from within Julia.
The package documentation, which is built through continuous integration (CI), includes runnable tutorial examples that are both available online for consultation and as Jupyter notebooks.
Continuous integration through GitHub actions also includes a fairly complete suite of tests.







# Statement of need





In spite of reduced realism, simplified ocean circulation models are at the forefront of oceanographic research, as evidenced by their increasing use in recent literature.
The low computational costs of steady-state circulation models allow for efficient optimization and inference/estimation of biogeochemical parameters, e.g., for the carbon cycle [@DeVries:2014; @Weber_etal_PNAS_2016; @Roshan_DeVries_NatCom_2017; @DeVries_etal_Nature_2017], for other macronutrient cycles (e.g. phosphorus, nitrogen, silicon) [@Kwon_Primeau_2006; @Kwon_Primeau_2008; @DeVries_etal_GBC_2013; @Holzer_etal_JGRO_2014; @Holzer_Brzezinski_GBC_2015; @DeVries_Weber_GBC_2017; @Wang_etal_Nature_2019], for stoichiometric ratios [@DeVries_Deutsch_NatGeo_2014; @Teng_etal_NatGeo_2014; @Weber_Deutsch_Nature_2010], and for micronutrients and trace metals such as dissolved iron [@Frants_etal_2016; @Pasquier_Holzer_2017], zinc [@Roshan_etal_GBC_2018], and argon [@Holzer_etal_GRL_2019], and hydrothermal gases such as helium [@Holzer_etal_EPSL_2017].
Efficient computation also allows for a wide range of experiments [@Primeau_etal:2013; @DeVries:2014; @Holzer_etal_GBC_2019].
Further, these models are amenable to the development of new metrics, thanks to their matrix representation, which enables the use of powerful linear-algebra techniques to efficiently compute exact partitions, pathways, timescales, eigen modes, moments, and other novel diagnostics [@Primeau_JPO_2005; @Holzer_Primeau_JGRO_2013; @Pasquier_Holzer_2016; @Holzer_etal_2016; @Fu_etal_2018; @Pasquier_Holzer_BG_2018; @Holzer_etal_JGRO_2020; @Holzer_etal_GBC_2021].
Finally, steady-state (and cyclo-stationary) circulation models are useful to bypass the need for long and costly spinups of large ocean general circulation models [@Khatiwala_etal_OM_2005; @Khatiwala_GBC_2007; @Khatiwala_OM_2008; @Khatiwala_Nature_2009; @Li_Primeau_OM_2008; @Bardin_etal_OM_2014; @Bardin_etal_OM_2016; @Huang_etal_2021].



Most of these studies relied on private implementations and required significant expertise.
Most often than not, these studies are hard to reproduce, if possible at all.
Comparisons between different circulation models are moreover complicated by the lack of standardization across models.
There is a need to facilitate and standardize the use of steady-state ocean-circulation models by providing
(i) an integrated tool for handling a number of different ocean-circulation models with tools for swapping circulations (including interpolating from one model grid to another),
(ii) a user-friendly interface for translating conceptual models of biogeochemical cycles into  the corresponding code (e.g., for sources, sinks, and vertical transport of tracers), and
(iii) solvers for efficient simulations, optimization, diagnosis, and statistical analysis.


To the best of our knowledge, AIBECS.jl is, at present, the best candidate for answering this need.
Although it serves a similar purpose, the fairly recently developed AWESOME OCIM [@John_etal_2020] currently lacks important features.
It provides only a single circulation model (the OCIM v1.0), cannot handle multiple tracers or nonlinear systems out of the box, and requires a MATLAB license.
In comparison, AIBECS.jl provides a free, open-source, unified framework for biogeochemical-tracer-modelling studies that use steady-state circulation models.
Among other advantages, AIBECS.jl offers more features, better computational efficiency, enhanced versatility, composability with other Julia packages, and ease of reproducibility (granted by version control and Julia's package manager) and improved syntax, which are pillars of modern scientific dissemination.
AIBECS.jl has been used for teaching and is currently used for research by members of the Marine Trace Element Laboratory of the University of Southern California.



The publication of circulation models as transport matrices from existing general circulation models [@Khatiwala_etal_OM_2005; @Khatiwala_GBC_2007; @Bardin_etal_OM_2014; @Bardin_etal_OM_2016; @Kvale_etal_GMD_2017; @Zanna_etal_PNAS_2019], and hopefully, future publication of transport matrix estimates from standard ocean general circulation models [e.g., @Chamberlain_etal_OM_2019] will increase the collection of circulations available from AIBECS.jl.
(Of particular interest to the broader community would be transport matrix models from the Climate Model Intercomparison Project (CMIP), which includes past and future simulations of the ocean circulation.)


Time and funding permitting, the scope of AIBECS could be extended to non-steady simulation through better integration with external state-of-the-art solvers, optimizers, and estimation algorithms provided, e.g.,  by the SciML ecosystem [@Rackauckas_Nie_JORS_2017].
This might allow for the eventual addition of cyclo-stationary circulations [e.g., CYCLOCIM: @Huang_etal_2021], which require adequate external Newton-Krylov solvers not yet available through AIBECS.jl.
Finally, bridging packages could be implemented for improved composability with, e.g., statistical ([Turing.jl](https://github.com/TuringLang/Turing.jl) [@ge2018t]), optimization, and plotting software ([Makie.jl](https://github.com/JuliaPlots/Makie.jl) [@Makie.jl]) to meet demand.

# Acknowledgements

FP and BP acknowledge funding from the Department of Energy grant DE-SC0016539 and from the National Science Foundation grant 1658380.
BP and SJ acknowledge funding from XXX.

# References