using DataDeps: standard_loadpath
using MD5: md5

"""
    _invalidate_stale_cache(name, filename, expected_md5)

Internal helper that closes a gap in DataDeps' caching: once a dep is
resolved on disk, `@datadep_str` short-circuits without re-checking the
hash, so bumping a registered URL/MD5 in source has no effect on machines
with a stale cache. This walks the DataDeps load path and, if a cached
`name/filename` exists with an MD5 that no longer matches `expected_md5`,
removes the directory so the next resolve triggers a fresh download.

Single-file deps only — tarball/multi-file deps need a different strategy.
"""
function _invalidate_stale_cache(
        name::AbstractString,
        filename::AbstractString,
        expected_md5::AbstractString;
        legacy_warning::Union{Nothing, AbstractString} = nothing,
    )
    for base in standard_loadpath
        cache_dir = joinpath(base, name)
        cache_file = joinpath(cache_dir, filename)
        if isfile(cache_file)
            actual = bytes2hex(open(md5, cache_file))
            if actual != expected_md5
                @warn "Cached file at $cache_file is stale (md5 $actual ≠ registered $expected_md5). Removing $cache_dir to trigger a fresh download."
                rm(cache_dir; recursive = true, force = true)
            end
            return
        elseif isdir(cache_dir) && legacy_warning !== nothing
            # Cache populated by a pre-0.17 unpack that did not preserve the
            # archive — we have nothing to hash, so we can only warn.
            @warn replace(legacy_warning, "{cache_dir}" => cache_dir)
            return
        end
    end
    return
end
