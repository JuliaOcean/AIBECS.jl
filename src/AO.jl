module AO

using DataDeps              # For storage location of data

function fallback_download(remotepath, localdir)
    @assert(isdir(localdir))
    filename = basename(remotepath)  # only works for URLs with filename as last part of name
    localpath = joinpath(localdir, filename)
    Base.download(remotepath, localpath)
    return localpath
end


function register_AO_files()
    register(
        DataDep(
            "AWESOME-OCIM",
            """
            References for OCIM1:
            - John, S. G., Liang, H., Weber, T., Devries, T., Primeau, F., Moore, K., Holzer, M., Mahowald, N., Gardner, W., Mishonov, A., Richardson, M., J., Faugere, Y., and Taburet, G. (2020). AWESOME OCIM: A simple, flexible, and powerful tool for modeling elemental cycling in the oceans. Chemical Geology, 533, 119403. doi: 10.1016/j.chemgeo.2019.119403.
            """,
            "https://github.com/hengdiliang/AWESOME-OCIM/archive/refs/heads/master.zip",
            sha2_256,
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
function download_and_unpack()
    register_AO_files()
    AO_path = datadep"AWESOME-OCIM"
    @info """You are about to download (and unpack) the AWESOME OCIM (AO)
          files as a zip file from GitHub and unpack them into:

              $AO_path

          You can run

              tree $AO_path

          (in shell mode) to check its contents.

          Please check with Seth John (sethjohn@usc.edu) for references 
          to cite if you use the other data contained in the AO files,
          e.g., GEOTRACES, WOA, Weber and John, and so on.

          Also note that downloading files like this from GitHub is 
          probably not very robust. However, some of the data files in 
          the AO repository do not exist anywhere else officially,
          therefore this seems like the best solution at this stage. 
          PRs welcome to improve this!
          """
    return nothing
end

end # module

export AO