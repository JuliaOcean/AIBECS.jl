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
using StringDistances
@reexport using OceanGrids            # To store the grid

function fallback_download(remotepath, localdir)
    @assert(isdir(localdir))
    filename = basename(remotepath)  # only works for URLs with filename as last part of name
    localpath = joinpath(localdir, filename)
    Downloads.download(remotepath, localpath)
    return localpath
end


# Fuzzy name matching
const DATASET_NAMES = ("Chien", "Kok")
const DATASET_LONGNAMES = Dict(
    "Chien" => "Chien_etal_2016",
    "Kok"   => "Kok_etal_2021"
)
rename(s) = findnearest(s, DATASET_NAMES, Levenshtein())[1]
renamelong(s) = DATASET_LONGNAMES[rename(s)]

# Aeolian source URLs
const URLs = Dict(
    "Chien" => "http://www.geo.cornell.edu/eas/PeoplePlaces/Faculty/mahowald/dust/Chienetal2016/post.aerosols.2x2.seasonal.nc",
    "Kok" => "https://dustcomm.atmos.ucla.edu/data/K21b/DustCOMM_source_region_wetdep_annual_PM20_abs.nc"
)
url(s) = URLs[rename(s)]

# Hashes
const SHAs = Dict(
    "Chien" => "33b0036bcca87019e0636ec2d59e910be989693819a2d57f2e4ab053a564101f", # On  2 Jun 2020
    "Kok"   => "65fec31419f75064a22a2c65247afc48e1f72218753605b463447a9933c378aa"  # On 14 Jul 2021
)
sha(s) = SHAs[rename(s)]

datadepname(s) = string("AIBECS-", renamelong(s))

# Create registry entry
function register_AeolianSource(s="Chien")
    register(
        DataDep(
            datadepname(s),
            """
            References:

            $(citations(s))
            """,
            url(s),
            sha(s),
            fetch_method = fallback_download
        )
    )
    return nothing
end

"""
    load()

Returns the 2D aeolian deposition fields from Chien et al. (2016).

You can specify a different dataset via, e.g.,

```julia
load("Kok")
```

At this stage, only two datasets are available:
- `"Chien"` (default) for different dust types (fires, biofuels, dust, ...)
- `"Kok"` for dust from different regions of origin
"""
load(s="Chien") = rename(s) == "Chien" ? load_Chien() : load_Kok()



#=============================#
#     Chien et al data        #
#=============================#
const Chien_AEROSOLTYPE_NAMES = (
    :fire,    # Fires
    :bioF,    # Biofuels
    :FF,      # FossilFuel
    :dust,    # Dust
    :salt,    # Seasalts
    :bioP,    # PrimaryBiogenicParticles
    :volc,    # Volcanoes
)
function load_Chien()
    s = "Chien"
    register_AeolianSource(s)
    nc_file = @datadep_str joinpath(datadepname(s), "post.aerosols.2x2.seasonal.nc")
    s_A_2D = Dataset(nc_file, "r") do ds
        Dict(
            :lat => ds["lat"][:],
            :lon => ds["lon"][:],
            (t => ds["dep"][:,:,:,i] for (i,t) in enumerate(Chien_AEROSOLTYPE_NAMES))...
        )
    end # ds is closed
    @info """You are about to use the Chien et al. (2016) data for aeolian deposition.
          If you use it for research, please cite:

          $(citations(s))

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the $(renamelong(s)) key.)
          """
    return s_A_2D
end



#=============================#
#       Kok et al data        #
#=============================#
const Kok_REGIONS_NAMES = (
    :NWAf,    # western North Africa
    :NEAf,    # eastern North Africa
    :Sahel,   # southern Sahara and Sahel
    :MECA,    # Middle East and central Asia
    :EAsia,   # East Asia
    :NAm,     # North America
    :Aus,     # Australia
    :SAm,     # South America
    :SAf,     # southern Africa
)

function load_Kok()
    s = "Kok"
    register_AeolianSource(s)
    nc_file = @datadep_str joinpath(datadepname(s), "DustCOMM_source_region_wetdep_annual_PM20_abs.nc")
# Units                = Total deposition (kg/m2/year)
    s_A_2D = Dataset(nc_file, "r") do ds
        Dict(
            :lat => ds["lat"][:],
            :lon => ds["lon"][:],
            (r => ds["Mean"][i,:,:] for (i,r) in enumerate(Kok_REGIONS_NAMES))...
        )
    end # ds is closed
    @info """You are about to use the Kok et al. (2021) data for aeolian deposition.
          If you use it for research, please cite:

          $(citations(s))

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the $(renamelong(s)) key.)
          """
    return s_A_2D
end



#=============================#
#          Citations          #
#=============================#
const CITATIONS = Dict(
    "Chien" => "- Chien, C.-T., K. R. M. Mackey, S. Dutkiewicz, N. M. Mahowald, J. M. Prospero, and A. Paytan (2016), Effects of African dust deposition on phytoplankton in the western tropical Atlantic Ocean off Barbados, Global Biogeochem. Cycles, 30, doi:10.1002/2015GB005334.",
    "Kok" => """
    - Description: Kok, J. F., Adebiyi, A. A., Albani, S., Balkanski, Y., Checa-Garcia, R., Chin, M., Colarco, P. R., Hamilton, D. S., Huang, Y., Ito, A., Klose, M., Li, L., Mahowald, N. M., Miller, R. L., Obiso, V., Pérez García-Pando, C., Rocha-Lima, A., and Wan, J. S.: Contribution of the world's main dust source regions to the global cycle of desert dust, Atmos. Chem. Phys., 21, 8169–8193, https://doi.org/10.5194/acp-21-8169-2021, 2021.
    - Dataset: Adebiyi, A. A., Kok, J. F., Wang, Y., Ito, A., Ridley, D. A., Nabat, P., and Zhao, C.: Dust Constraints from joint Observational-Modelling-experiMental analysis (DustCOMM): comparison with measurements and model simulations, Atmos. Chem. Phys., 20, 829–863, https://doi.org/10.5194/acp-20-829-2020, 2020.
    - Source apportionment framework: Kok, J. F., Adebiyi, A. A., Albani, S., Balkanski, Y., Checa-Garcia, R., Chin, M., Colarco, P. R., Hamilton, D. S., Huang, Y., Ito, A., Klose, M., Leung, D. M., Li, L., Mahowald, N. M., Miller, R. L., Obiso, V., Pérez García-Pando, C., Rocha-Lima, A., Wan, J. S., and Whicker, C. A.: Improved representation of the global dust cycle using observational constraints on dust properties and abundance, Atmos. Chem. Phys., 21, 8127–8167, https://doi.org/10.5194/acp-21-8127-2021, 2021.
    - Climate model used: Scanza, R. A., Hamilton, D. S., Perez Garcia-Pando, C., Buck, C., Baker, A., and Mahowald, N. M.: Atmospheric processing of iron in mineral and combustion aerosols: development of an intermediate-complexity mechanism suitable for Earth system models, Atmos. Chem. Phys., 18, 14175–14196, https://doi.org/10.5194/acp-18-14175-2018, 2018.
    """
)
citations(s="Chien") = CITATIONS[rename(s)]

end # end module

export AeolianSources


