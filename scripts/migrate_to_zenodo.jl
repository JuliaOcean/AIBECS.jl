#!/usr/bin/env julia
# Migrate the AIBECS OCIM/OCCA datasets from FigShare to Zenodo.
#
# Run:
#     export ZENODO_TOKEN=...      # personal access token, scope: deposit:write
#     julia --project=scripts scripts/migrate_to_zenodo.jl
#
# What it does, per family (OCIM0, OCIM1, OCCA, OCIM2_48L, OCIM2):
#   1. Downloads each file from FigShare into a temp dir.
#   2. Verifies the FigShare MD5 against the value baked into the AIBECS source
#      (so we know we are uploading the exact bytes users already trust).
#   3. Creates a Zenodo *draft* deposition for the family with citation metadata.
#   4. Uploads every file to the deposition's bucket URL.
#   5. Reads the resulting per-file MD5 and download URL from the Zenodo API.
#
# Drafts are NOT published. Open each one in the Zenodo web UI to review and
# click Publish. The script can be re-run safely: existing drafts are reused
# (matched by family name in the deposition title) and files already present
# in the bucket are skipped.
#
# At the end the script prints:
#   - a Julia snippet ready to paste into src/OCIM*.jl
#   - a markdown table fragment ready to drop into docs/src/reference/datasets.md

using Downloads
using HTTP
using JSON3
using MD5
using Printf

const ZENODO = "https://zenodo.org/api"
const TOKEN = get(ENV, "ZENODO_TOKEN", "")
isempty(TOKEN) && error("Set ZENODO_TOKEN before running this script.")

const AUTH_HEADER = ["Authorization" => "Bearer $TOKEN"]
const JSON_HEADERS = ["Content-Type" => "application/json"; AUTH_HEADER]

# ---------------------------------------------------------------------------
# Source of truth: mirrored from src/OCIM*.jl. Update here if upstream changes.
# ---------------------------------------------------------------------------

const CREATORS = [Dict(
    "name"        => "Pasquier, Benoît",
    "affiliation" => "School of Mathematics and Statistics, University of New South Wales",
    "orcid"       => "0000-0002-3838-5976",
)]

# Citation strings copied verbatim from the existing AIBECS modules.
const CITATIONS = Dict(
    "OCIM0" => """
    - Primeau, F. W., Holzer, M., and DeVries, T. (2013), Southern Ocean nutrient trapping and the efficiency of the biological pump, J. Geophys. Res. Oceans, 118, 2547–2564, doi:10.1002/jgrc.20181.
    - DeVries, T. and F. Primeau, 2011: Dynamically and Observationally Constrained Estimates of Water-Mass Distributions and Ages in the Global Ocean. J. Phys. Oceanogr., 41, 2381–2401, doi:10.1175/JPO-D-10-05011.1
    """,
    "OCIM1" => """
    - DeVries, T., 2014: The oceanic anthropogenic CO2 sink: Storage, air-sea fluxes, and transports over the industrial era, Global Biogeochem. Cycles, 28, 631–647, doi:10.1002/2013GB004739.
    - DeVries, T. and F. Primeau, 2011: Dynamically and Observationally Constrained Estimates of Water-Mass Distributions and Ages in the Global Ocean. J. Phys. Oceanogr., 41, 2381–2401, doi:10.1175/JPO-D-10-05011.1
    """,
    "OCCA" => "Forget, G., 2010: Mapping Ocean Observations in a Dynamical Framework: A 2004–06 Ocean Atlas. J. Phys. Oceanogr., 40, 1201–1221, https://doi.org/10.1175/2009JPO4043.1",
    "OCIM2" => "DeVries, T., & Holzer, M. (2019). Radiocarbon and helium isotope constraints on deep ocean ventilation and mantle-3He sources. Journal of Geophysical Research: Oceans, 124, 3036–3057. https://doi.org/10.1029/2018JC014716",
    "OCIM2_48L" => """
    - Holzer, M., DeVries, T. & de Lavergne, C. Diffusion controls the ventilation of a Pacific Shadow Zone above abyssal overturning. Nat Commun 12, 4348 (2021). https://doi.org/10.1038/s41467-021-24648-x
    - DeVries, T., & Holzer, M. (2019). Radiocarbon and helium isotope constraints on deep ocean ventilation and mantle-3He sources. Journal of Geophysical Research: Oceans, 124, 3036–3057. https://doi.org/10.1029/2018JC014716
    """,
)

const TITLES = Dict(
    "OCIM0"     => "OCIM0 ocean circulation matrix and grid (JLD2)",
    "OCIM1"     => "OCIM1 ocean circulation matrix and grid (JLD2)",
    "OCCA"      => "OCCA ocean circulation matrix and grid (JLD2)",
    "OCIM2"     => "OCIM2 ocean circulation matrices and grids (JLD2, 11 variants)",
    "OCIM2_48L" => "OCIM2-48L ocean circulation matrix and grid (tar.gz)",
)

const DESCRIPTIONS = Dict(
    "OCIM0" => """
    JLD2 conversion of the OCIM v0 transport matrix and grid for use with
    AIBECS.jl. Built from the upstream MATLAB distribution by the script
    in <https://github.com/briochemc/OceanCirculations>. Migrated from
    FigShare to Zenodo for hosting reliability.
    """,
    "OCIM1" => """
    JLD2 conversion of the OCIM1 (CTL) transport matrix and grid for use
    with AIBECS.jl. Built from the upstream MATLAB distribution by the
    script in <https://github.com/briochemc/OceanCirculations>. Migrated
    from FigShare to Zenodo for hosting reliability.
    """,
    "OCCA" => """
    JLD2 conversion of the OCCA transport matrix and grid for use with
    AIBECS.jl. Built from the upstream MATLAB distribution by the script
    in <https://github.com/briochemc/OceanCirculations>. Migrated from
    FigShare to Zenodo for hosting reliability.
    """,
    "OCIM2" => """
    JLD2 conversions of the OCIM2 transport matrices and grids (11 variants
    spanning Ki, Kv, and helium settings) for use with AIBECS.jl. Built
    from the upstream MATLAB distribution by the script in
    <https://github.com/briochemc/OceanCirculations>. Migrated from
    FigShare to Zenodo for hosting reliability.
    """,
    "OCIM2_48L" => """
    Tar.gz bundle of the OCIM2-48L (48-layer) transport matrix and grid
    for use with AIBECS.jl. Migrated from FigShare to Zenodo for hosting
    reliability.
    """,
)

# Each entry: (filename, ndownloader URL, src/ MD5, FigShare article ID)
# Article IDs come from briochemc's FigShare account; one article per JLD2 file.
const FAMILIES = [
    ("OCIM0", [
        ("OCIM0.1.jld2", "https://ndownloader.figshare.com/files/28336086", "458a463a9f1fdc223a6cfb025daa2d47", 8317085),
    ]),
    ("OCIM1", [
        ("OCIM1_CTL.jld2", "https://ndownloader.figshare.com/files/28335882", "eaa57b42e7edec0fe965575e9938c66d", 8241176),
    ]),
    ("OCCA", [
        ("OCCA.jld2", "https://ndownloader.figshare.com/files/28336173", "f2a2a2d295f85771c2302e6c0eb35f4c", 12386498),
    ]),
    ("OCIM2_48L", [
        ("OCIM2_48L_base.tar.gz", "https://ndownloader.figshare.com/files/28468077", "7938f20f06eff072b2791b571d2bb9d7", 14802732),
    ]),
    ("OCIM2", [
        ("OCIM2_CTL_He.jld2",              "https://ndownloader.figshare.com/files/28336284", "da15192381ef04e9f7cce7c886eb4833", 11911104),
        ("OCIM2_CTL_noHe.jld2",            "https://ndownloader.figshare.com/files/28336299", "abede7e16b75cb3f9b7d25367772455e", 11911122),
        ("OCIM2_KiHIGH_He.jld2",           "https://ndownloader.figshare.com/files/28336302", "3c8312f9aeb9bb6cb6ac32a8d84fb901", 11911125),
        ("OCIM2_KiHIGH_noHe.jld2",         "https://ndownloader.figshare.com/files/28336311", "781b8c058eb74a81d10befa83faf0d08", 11911128),
        ("OCIM2_KiLOW_He.jld2",            "https://ndownloader.figshare.com/files/28336317", "ee6c6f14c6bd90c1980daf8a220f0551", 11911131),
        ("OCIM2_KiLOW_noHe.jld2",          "https://ndownloader.figshare.com/files/28336323", "47eabd73588739344c81f72d39367f27", 11911134),
        ("OCIM2_KvHIGH_He.jld2",           "https://ndownloader.figshare.com/files/28336326", "ecd18348d062fb1c1e1518244b2b8275", 11911137),
        ("OCIM2_KvHIGH_KiHIGH_noHe.jld2",  "https://ndownloader.figshare.com/files/28336329", "c71514b13802c0e497f033612f586f5a", 11911140),
        ("OCIM2_KvHIGH_KiLOW_He.jld2",     "https://ndownloader.figshare.com/files/28336332", "99edb2894a2ee274e68f3309da1ed200", 11911143),
        ("OCIM2_KvHIGH_KiLOW_noHe.jld2",   "https://ndownloader.figshare.com/files/28336341", "ff8fb1a0add22bde6f3fed393fa7336f", 11911146),
        ("OCIM2_KvHIGH_noHe.jld2",         "https://ndownloader.figshare.com/files/28336353", "0c0402b3796c3a7d14c0bdc2db7f0aa7", 11911149),
    ]),
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

hexmd5(path) = bytes2hex(open(io -> md5(io), path))

# Files where the FigShare-API MD5 differs from what's currently in src/OCIM*.jl.
# Key = filename, Value = (src_md5, figshare_api_md5). Reported at end.
const STALE_HASHES = Dict{String, Tuple{String, String}}()
record_stale!(name, src_md5, fs_md5) = (STALE_HASHES[name] = (src_md5, fs_md5))

# Ask FigShare itself what MD5 it has for `filename` inside `article_id`.
# Returns the MD5 string, or nothing if the article/file can't be resolved.
function figshare_api_md5(article_id, filename)
    try
        r = HTTP.get("https://api.figshare.com/v2/articles/$article_id"; status_exception=false)
        r.status == 200 || return nothing
        body = JSON3.read(r.body)
        for f in get(body, :files, [])
            string(f.name) == filename && return string(f.computed_md5)
        end
        @warn "  file not found in FigShare article" article_id filename
        return nothing
    catch e
        @warn "  FigShare API lookup failed" article_id exception=e
        return nothing
    end
end

function fetch_with_check(url, dest, src_md5, article_id; max_attempts=4)
    name = basename(dest)
    fs_md5 = figshare_api_md5(article_id, name)

    if fs_md5 === nothing
        @warn "  FigShare API didn't return an MD5; falling back to src/ MD5" file=name
    elseif fs_md5 != src_md5
        @warn "  src/ MD5 is stale (FigShare's canonical MD5 differs)" file=name src_md5 figshare_md5=fs_md5
        record_stale!(name, src_md5, fs_md5)
    end
    target_md5 = something(fs_md5, src_md5)  # what we'll actually verify against

    if isfile(dest) && filesize(dest) > 0
        cached = hexmd5(dest)
        if cached == target_md5
            @info "  cached, verified" file=name md5=cached size=filesize(dest)
            return
        end
        @info "  cached but MD5 mismatch, will redownload" file=name cached target=target_md5
    end

    seen = String[]
    for attempt in 1:max_attempts
        @info "  downloading from FigShare" url attempt
        try
            Downloads.download(url, dest)
        catch e
            @warn "  download error" attempt exception=e
            attempt == max_attempts && rethrow()
            sleep(2.0 * attempt)
            continue
        end
        sz = filesize(dest)
        got = hexmd5(dest)
        push!(seen, got)
        if got == target_md5
            @info "  download verified against $(fs_md5 === nothing ? "src/" : "FigShare-API") MD5" file=name size=sz
            return
        end
        @warn "  download MD5 differs from authoritative target, retrying" attempt target=target_md5 got=got size=sz
        sleep(1.0 * attempt)
    end
    error("MD5 mismatch for $dest after $max_attempts attempts: target=$target_md5 observed=$(seen)")
end

function zenodo_get(path)
    r = HTTP.get(ZENODO * path; headers=AUTH_HEADER)
    return JSON3.read(r.body)
end

function zenodo_post(path, body)
    r = HTTP.post(ZENODO * path; headers=JSON_HEADERS, body=JSON3.write(body))
    return JSON3.read(r.body)
end

function zenodo_put_json(path, body)
    r = HTTP.put(ZENODO * path; headers=JSON_HEADERS, body=JSON3.write(body))
    return JSON3.read(r.body)
end

function find_existing_draft(title)
    # Search the user's depositions for a draft with this title.
    r = HTTP.get(ZENODO * "/deposit/depositions?status=draft&size=100"; headers=AUTH_HEADER)
    for d in JSON3.read(r.body)
        string(get(d, :title, "")) == title && return d
    end
    return nothing
end

function ensure_draft(family)
    title = TITLES[family]
    existing = find_existing_draft(title)
    if existing !== nothing
        @info "  reusing existing draft" id=existing.id title
        # The list endpoint returns abbreviated records, so re-fetch to get bucket URL.
        return zenodo_get("/deposit/depositions/$(existing.id)")
    end
    @info "  creating new draft" title
    body = Dict("metadata" => Dict(
        "upload_type" => "dataset",
        "title"       => title,
        "creators"    => CREATORS,
        "description" => strip(DESCRIPTIONS[family]),
        "notes"       => strip(CITATIONS[family]),
        "keywords"    => ["AIBECS.jl", "ocean circulation", "transport matrix", "biogeochemistry"],
        "access_right"=> "open",
        "license"     => "cc-by-4.0",
    ))
    return zenodo_post("/deposit/depositions", body)
end

function upload_file(bucket_url, filepath)
    name = basename(filepath)
    bytes = read(filepath)
    HTTP.put("$bucket_url/$name"; headers=AUTH_HEADER, body=bytes)
    @info "  uploaded" name length=length(bytes)
end

function existing_filenames(deposition_id)
    files = zenodo_get("/deposit/depositions/$deposition_id/files")
    return Set(string(f.filename) for f in files)
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function migrate()
    workdir = mktempdir(; prefix="aibecs_migration_")
    @info "Working directory" workdir

    results = Dict{String, Any}()  # family => Vector of (filename, zenodo_url, md5)

    for (family, files) in FAMILIES
        @info "▶ Family: $family"
        family_dir = joinpath(workdir, family)
        mkpath(family_dir)

        # 1. Download + verify (against FigShare API MD5; falls back to src/ MD5)
        for (name, url, src_md5, article_id) in files
            fetch_with_check(url, joinpath(family_dir, name), src_md5, article_id)
        end

        # 2. Ensure draft deposition exists
        draft = ensure_draft(family)
        bucket = string(draft.links.bucket)
        already = existing_filenames(draft.id)

        # 3. Upload missing files
        for f in files
            name = f[1]
            if name in already
                @info "  already present, skipping upload" name
                continue
            end
            upload_file(bucket, joinpath(family_dir, name))
        end

        # 4. Read back per-file URL + MD5 from API
        rows = []
        for f in zenodo_get("/deposit/depositions/$(draft.id)/files")
            push!(rows, (
                name = string(f.filename),
                url  = "https://zenodo.org/records/$(draft.id)/files/$(f.filename)",
                md5  = string(f.checksum),
                size = Int(f.filesize),
            ))
        end
        results[family] = (deposition_id = draft.id, html_url = string(draft.links.html), rows = rows)
        @info "✔ $family draft ready" review_url=string(draft.links.html)
    end

    return results
end

# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

function variant_from_filename(family, name)
    family == "OCIM2" || return nothing
    base = replace(name, "OCIM2_" => "")
    base = replace(base, ".jld2" => "")
    return base
end

function print_julia_snippets(results)
    println("\n# ===== Paste into src/OCIM*.jl =====\n")
    for family in ("OCIM0", "OCIM1", "OCCA", "OCIM2_48L", "OCIM2")
        haskey(results, family) || continue
        rows = results[family].rows
        println("# $family")
        if family == "OCIM2"
            println("const OCIM2_URLs = Dict(")
            for r in rows
                v = variant_from_filename(family, r.name)
                @printf("    %-26s => \"%s\",\n", "\"$v\"", r.url)
            end
            println(")")
            println("const OCIM2_MD5 = Dict(")
            for r in rows
                v = variant_from_filename(family, r.name)
                @printf("    %-26s => \"%s\",\n", "\"$v\"", r.md5)
            end
            println(")")
        else
            r = rows[1]
            const_name = family == "OCIM2_48L" ? "URL" : "$(family)_URL"
            md5_const  = family == "OCIM2_48L" ? "OCIM2_48L_MD5" : "$(family)_MD5"
            println("const $const_name = \"$(r.url)\"")
            println("const $md5_const = \"$(r.md5)\"")
        end
        println()
    end
end

function print_markdown_table_fragment(results)
    println("\n# ===== Paste into docs/src/reference/datasets.md =====\n")
    println("| Module | Variant | Grid size | Citation | Source | Size (MB) |")
    println("|---|---|---|---|---|---|")
    grid = Dict(
        "OCIM0" => "180×90×24", "OCIM1" => "180×91×24", "OCCA" => "180×80×10",
        "OCIM2" => "180×91×24", "OCIM2_48L" => "180×91×48",
    )
    cite = Dict(
        "OCIM0"     => "[DeVries & Primeau (2011)](https://doi.org/10.1175/JPO-D-10-05011.1); [Primeau et al. (2013)](https://doi.org/10.1002/jgrc.20181)",
        "OCIM1"     => "[DeVries (2014)](https://doi.org/10.1002/2013GB004739)",
        "OCCA"      => "[Forget (2010)](https://doi.org/10.1175/2009JPO4043.1)",
        "OCIM2"     => "[DeVries & Holzer (2019)](https://doi.org/10.1029/2018JC014716)",
        "OCIM2_48L" => "[Holzer, DeVries & de Lavergne (2021)](https://doi.org/10.1038/s41467-021-24648-x)",
    )
    for family in ("OCIM0", "OCIM1", "OCIM2", "OCIM2_48L", "OCCA")
        haskey(results, family) || continue
        for r in results[family].rows
            variant = if family == "OCIM2"
                v = variant_from_filename(family, r.name)
                v == "CTL_He" ? "CTL_He (default)" : v
            elseif family == "OCIM1"
                "CTL"
            elseif family == "OCIM2_48L"
                "base"
            else
                "(single)"
            end
            mb = round(Int, r.size / 1e6)
            println("| `$family` | $variant | $(grid[family]) | $(cite[family]) | [link]($(r.url)) | $mb |")
        end
    end
end

# ---------------------------------------------------------------------------

results = migrate()
print_julia_snippets(results)
print_markdown_table_fragment(results)

println("\nReview each draft and click Publish:")
for (family, info) in results
    println("  $family => $(info.html_url)")
end

if !isempty(STALE_HASHES)
    println("\n# ===== Stale MD5s detected in src/ =====")
    println("# These files on FigShare have a different MD5 than what's hardcoded")
    println("# in the AIBECS source. The Zenodo upload uses FigShare's actual bytes,")
    println("# so the new (md5, hash) tuples printed above already reflect the correct")
    println("# values. The src/ hashes were likely updated incorrectly at some point.")
    for (name, (expected, actual)) in sort(collect(STALE_HASHES))
        println("  $name: src=$expected  figshare=$actual")
    end
end
