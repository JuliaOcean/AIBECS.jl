using AIBECS
using Unitful

# ---------------------------------------------------------------------------
# Circulation tiers. Each entry maps a tier to a list of (label, loader)
# tuples. The loader returns (grd, T_circ) — the grid and the transport
# matrix with units stripped.
# ---------------------------------------------------------------------------

_strip(load) = function ()
    grd, T_unit = load()
    return grd, ustrip.(T_unit)
end

const TIERS = (
    small  = [("OCIM0", _strip(OCIM0.load)),
              ("OCCA",  _strip(OCCA.load))],
    medium = [("OCIM1", _strip(OCIM1.load)),
              ("OCIM2", _strip(OCIM2.load))],
    large  = [("OCIM2_48L", _strip(OCIM2_48L.load))],
)
