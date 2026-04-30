#---------------------------------------------------------
# # [AWESOME OCIM toolbox](@id awesome-ocim)
#---------------------------------------------------------

# The [AWESOME OCIM (AO)](https://github.com/hengdiliang/AWESOME-OCIM)
# is a MATLAB toolbox by John, Liang, Weber, DeVries, Primeau, Moore,
# Holzer, and Mahowald (2020) that bundles OCIM1 transport matrices
# alongside auxiliary GEOTRACES, WOA, and Weber-and-John datasets. AIBECS
# does not expose the AO contents directly, but it can fetch and unpack
# the toolbox locally so you can browse the underlying files yourself.

# ## Download and unpack

# `AO.download_and_unpack()` registers a `DataDep`, downloads the
# upstream GitHub zip on first call, and returns the path of the
# unpacked tree:

#md # ```julia
#md # using AIBECS
#md # AO_path = AO.download_and_unpack()
#md # ```

# We hide the actual call from this rendered page so the docs build does
# not keep re-downloading from GitHub on every CI run.

# ## What you get

# Inside `AO_path` you will find:
#
# - `OCIM1/` — MATLAB-format transport matrices and grid info,
# - `data/` — bundled observational fields (GEOTRACES, WOA, Weber–John),
# - `util/` — MATLAB helper scripts.

# Refer to the AO citation in the docstring of [`AO`](@ref) when using
# any of the bundled data.

# ## Caveat: GitHub-zip distribution

# The AO is served as a GitHub archive zip rather than from a stable
# DOI mirror, so downloads depend on GitHub's availability and the
# upstream branch staying intact. If you plan to use the AO data
# heavily, copy the unpacked tree into your own storage. PRs to point
# this at a more durable mirror are welcome.
