"""
This `OCIM2_48L` module is used to load the OCIM2-48L matrices and grid for use in AIBECS.

!!! tip
    To load the default OCIM2_48L matrix and grid, do
    ```
    julia> grd, T = OCIM2_48L.load()
    ```
    But you can also load the other matrices by specifying which version you want, e.g.,
    ```
    julia> grd, T = OCIM2_48L.load(version="OCIM2_48L_KiHIGH_noHe")
    ```
    See *DeVries and Holzer* (2019) for more details


!!! note
    The files, that are downloaded from a public and persistant URL in FigShare,
    were created with the code available at https://github.com/briochemc/OceanCirculations.
"""
module OCIM2_48L

using SparseArrays          # for sparse matrix
using DataDeps              # for automated download/storage of data
using Downloads
using MAT                   # for reading OCIM2_48L transport operator matrix from MAT file
using NCDatasets            # for reading OCIM2_48L grid from NetCDF file
using Unitful               # for units
using Unitful: s, yr, m, km, °
using Reexport
@reexport using OceanGrids            # To store the grid

function fallback_download(remotepath, localdir)
    @assert(isdir(localdir))
    filename = basename(remotepath)  # only works for URLs with filename as last part of name
    localpath = joinpath(localdir, filename)
    Downloads.download(remotepath, localpath)
    return localpath
end

# OCIM2_48L URL
url() = "https://files.figshare.com/28468077/OCIM2_48L_base.tar.gz"

# OCIM2_48L Hashes
sha() = "7c6e7981df69122957c9f2475be7ffe958819ff44d02be6d004886b6dfc3fac7"

OCIM2versionerror(version) = error("""`$version` is not a valid OCIM2 version name.

                                   Valid versions are `CTL_He`, `CTL_noHe`, `KiHIGH_He`, `KiHIGH_noHe`, `KiLOW_He`, `KiLOW_noHe`, `KvHIGH_He`, `KvHIGH_KiHIGH_noHe`, `KvHIGH_KiLOW_He`, `KvHIGH_KiLOW_noHe`, and `KvHIGH_noHe`.

                                   See *DeVries and Holzer* (2019) for more details on OCIM2 configurations.""")


# Create registry entry for OCIM in JLD2 format
function register_OCIM2_48L()
    register(
        DataDep(
            "AIBECS-OCIM2_48L",
            """
            References:
            - $(citations())
            """,
            url(),
            sha(),
            fetch_method = fallback_download,
            post_fetch_method = unpack
        )
    )
    return nothing
end

"""
    load

Returns the grid and the transport matrix of the OCIM2_48L model.

See *DeVries and Holzer* (2019) and *Holzer et al.* (2021) for more details.
"""
function load()
    register_OCIM2_48L()
    files_path = @datadep_str "AIBECS-OCIM2_48L"
    @info """You are about to use the OCIM2_48L model.
          If you use it for research, please cite:

          $(citations())

          You can find the corresponding BibTeX entries in the CITATION.bib file
          at the root of the AIBECS.jl package repository.
          (Look for the "DeVries_Holzer_2019" and "Holzer_etal_2021" keys.)
          """

    # Convert convergence T in yr⁻¹ from original MAT file to divergence T in s⁻¹
    T = -ustrip.(s^-1, matread(joinpath(files_path, "OCIM2_48L_base_transport.mat"))["TR"] * yr^-1);

    # Create grid object from NetCDF file variables
    # TODO: Add Hyrothermal He fluxes as for the OCIM2 load function
    grid = Dataset(joinpath(files_path, "OCIM2_48L_base_data.nc"), "r") do ds
        wet3D = ds["ocnmask"][:] .== 1
        nlat, nlon, ndepth = size(wet3D)
        lat_3D = ds["tlat"][:]°
        lon_3D = ds["tlon"][:]°
        depth_3D = ds["tz"][:]m
        lat = unique(lat_3D)
        lon = unique(lon_3D)
        depth = unique(depth_3D)
        ulat = unique(ds["ulon"][:])° # ulat↔ulon in OCIM2-48L original files
        ulon = unique(ds["ulat"][:])° # ulat↔ulon in OCIM2-48L original files
        depth_top_3D = ds["wz"][:]m
        depth_top = unique(depth_top_3D)
        δlat = 2(ulat - lat)
        δlon = 2(ulon - lon)
        δdepth = 2(depth - depth_top)
        R = 6371.0km
        δy = R * δlat ./ 360°
        δy_3D = repeat(reshape(δy, (nlat,1,1)), outer=(1,nlon,ndepth))
        A_3D = ds["area"][:]m^2
        δx_3D = A_3D ./ δy_3D
        volume_3D = ds["vol"]m^3
        δz_3D = volume_3D ./ A_3D
        A_2D = A_3D[:,:,1]
        nboxes = count(wet3D)
        OceanRectilinearGrid(
            lat,
            lon,
            depth,
            δlat,
            δlon,
            δdepth,
            lat_3D,
            lon_3D,
            depth_3D,
            δy,
            δx_3D,
            δy_3D,
            δz_3D,
            volume_3D,
            depth_top,
            depth_top_3D,
            A_2D,
            wet3D,
            nlon,
            nlat,
            ndepth,
            nboxes
        )
    end

    return grid, T

end

citations() = """
              - Holzer, M., DeVries, T. & de Lavergne, C. Diffusion controls the ventilation of a Pacific Shadow Zone above abyssal overturning. Nat Commun 12, 4348 (2021). https://doi.org/10.1038/s41467-021-24648-x
              - DeVries, T., & Holzer, M. (2019). Radiocarbon and helium isotope constraints on deep ocean ventilation and mantle-3He sources. Journal of Geophysical Research: Oceans, 124, 3036–3057. https://doi.org/10.1029/2018JC014716
              """

end # end module

export OCIM2_48L

