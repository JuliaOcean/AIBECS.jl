module AO

using DataDeps              # For storage location of data

function fallback_download(remotepath, localdir)
    @assert(isdir(localdir))
    filename = basename(remotepath)  # only works for URLs with filename as last part of name
    localpath = joinpath(localdir, filename)
    Base.download(remotepath, localpath)
    return localpath
end


function register_AO_files(AO_version)
    register(
        DataDep(
            string("AO_", AO_version, "_files"),
            """
            References for OCIM1:
            - DeVries, T., 2014: The oceanic anthropogenic CO2 sink: Storage, air‐sea fluxes, and transports over the industrial era, Global Biogeochem. Cycles, 28, 631–647, doi:10.1002/2013GB004739.
            - DeVries, T. and F. Primeau, 2011: Dynamically and Observationally Constrained Estimates of Water-Mass Distributions and Ages in the Global Ocean. J. Phys. Oceanogr., 41, 2381–2401, https://doi.org/10.1175/JPO-D-10-05011.1
            """,
            string("http://www.mtel.rocks/mtel/awesomeOCIM_files/AwesomeOCIM_", AO_version, ".zip"),
            "46b8bfad2b9c73e04cbb2688efd7b477af6ce61ad0340aa894c950c5ba53fa14",
            fetch_method = fallback_download,
            post_fetch_method = unpack
        )
    )
    return nothing
end

"""
    download_and_unpack

Downloads and unpacks the AO zip file from the MTEL website.
"""
function download_and_unpack(AO_version)
    register_AO_files(AO_version)
    @info """You are about to download (and unpack) the AWESOME OCIM (AO) $AO_version files as a zip file from the MTEL website (https://www.mtel.rocks).

          Please check with Seth John (sethjohn@usc.edu) for references to cite if you use the other data contained in the AO files, e.g., GEOTRACES, WOA, etc.
          (The references to DeVries et al. pertain to the OCIM v1 product, which is essentially the transport matrix representing the circulation.)

          """
    files_path = @datadep_str string("AO_", AO_version, "_files/")
    println("AO folder contains:\n",readdir(files_path)...)
    return nothing
end

end # module

export AO