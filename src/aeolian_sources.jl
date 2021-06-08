"""
This `AeolianSources` module is to create aeolian sources for use in AIBECS-generated models.
"""
module AeolianSources

using SparseArrays          # For sparse matrix
using DataDeps              # For storage location of data
using Downloads
using NCDatasets
using Unitful               # for units
using Reexport
using Interpolations
@reexport using OceanGrids            # To store the grid

function fallback_download(remotepath, localdir)
    @assert(isdir(localdir))
    filename = basename(remotepath)  # only works for URLs with filename as last part of name
    localpath = joinpath(localdir, filename)
    Downloads.download(remotepath, localpath)
    return localpath
end

# Aeolian source URLs
url() = "http://www.geo.cornell.edu/eas/PeoplePlaces/Faculty/mahowald/dust/Chienetal2016/post.aerosols.2x2.seasonal.nc"

# Hashes
sha() = "33b0036bcca87019e0636ec2d59e910be989693819a2d57f2e4ab053a564101f" # On 2 June 2020

# Create registry entry
function register_AeolianSource()
    register(
        DataDep(
            "AIBECS-Chienetal2016",
            """
            References:
            - $(citation())
            """,
            url(),
            sha(),
            fetch_method = fallback_download
        )
    )
    return nothing
end

"""
    load

Returns the 2D deposition fields.
"""
function load()
    register_AeolianSource()
    nc_file = @datadep_str string("AIBECS-Chienetal2016/", "post.aerosols.2x2.seasonal.nc")
    s_A_2D = Dataset(nc_file, "r") do ds
        Dict(
            :lat => ds["lat"][:],
            :lon => ds["lon"][:],
            :Fires                    => ds["dep"][:,:,:,1],
            :Biofuels                 => ds["dep"][:,:,:,2],
            :FossilFuel               => ds["dep"][:,:,:,3],
            :Dust                     => ds["dep"][:,:,:,4],
            :Seasalts                 => ds["dep"][:,:,:,5],
            :PrimaryBiogenicParticles => ds["dep"][:,:,:,6],
            :Volcanoes                => ds["dep"][:,:,:,7]
        )
    end # ds is closed
    @info """You are about to use the Chien et al. (2016) data for aeolian deposition.
          If you use it for research, please cite:

          - $(citation())

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "Chien_etal_2016" key.)
          """
    return s_A_2D
end

citation() = "Chien, C.-T., K. R. M. Mackey, S. Dutkiewicz, N. M. Mahowald, J. M. Prospero, and A. Paytan (2016), Effects of African dust deposition on phytoplankton in the western tropical Atlantic Ocean off Barbados, Global Biogeochem. Cycles, 30, doi:10.1002/2015GB005334."

end # end module

export AeolianSources


