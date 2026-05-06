module AIBECSOCIM2_48LExt

using AIBECS
using AIBECS.OCIM2_48L: OCIM2_48L, register_OCIM2_48L, CITATION
using DataDeps
using MAT
using NCDatasets
using OceanGrids
using Unitful
using Unitful: s, yr, m, km, °

function OCIM2_48L.load()
    register_OCIM2_48L()
    files_path = @datadep_str "AIBECS-OCIM2_48L"
    @info """You are about to use the OCIM2_48L model.
    If you use it for research, please cite:

    $CITATION

    You can find the corresponding BibTeX entries in the CITATION.bib file
    at the root of the AIBECS.jl package repository.
    (Look for the "DeVries_Holzer_2019" and "Holzer_etal_2021" keys.)
    """

    # Convert convergence T in yr⁻¹ from original MAT file to divergence T in s⁻¹
    T = -ustrip.(s^-1, matread(joinpath(files_path, "OCIM2_48L_base_transport.mat"))["TR"] * yr^-1)

    # Create grid object from NetCDF file variables
    # TODO: Add Hyrothermal He fluxes as for the OCIM2 load function
    grid = Dataset(joinpath(files_path, "OCIM2_48L_base_data.nc"), "r") do ds
        wet3D = Array(ds["ocnmask"]) .== 1
        nlat, nlon, ndepth = size(wet3D)
        lat_3D = Array(ds["tlat"])°
        lon_3D = Array(ds["tlon"])°
        depth_3D = Array(ds["tz"])m
        lat = unique(lat_3D)
        lon = unique(lon_3D)
        depth = unique(depth_3D)
        ulat = unique(Array(ds["ulon"]))° # ulat↔ulon in OCIM2-48L original files
        ulon = unique(Array(ds["ulat"]))° # ulat↔ulon in OCIM2-48L original files
        depth_top_3D = Array(ds["wz"])m
        depth_top = unique(depth_top_3D)
        δlat = 2(ulat - lat)
        δlon = 2(ulon - lon)
        δdepth = 2(depth - depth_top)
        R = 6371.0km
        δy = R * δlat ./ 360°
        δy_3D = repeat(reshape(δy, (nlat, 1, 1)), outer = (1, nlon, ndepth))
        A_3D = Array(ds["area"])m^2
        δx_3D = A_3D ./ δy_3D
        volume_3D = Array(ds["vol"])m^3
        δz_3D = volume_3D ./ A_3D
        A_2D = A_3D[:, :, 1]
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

end # module
