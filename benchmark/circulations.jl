using AIBECS
using Unitful

# ---------------------------------------------------------------------------
# Circulation tiers. Each entry maps a tier to a list of (label, loader)
# tuples. The loader returns (grd, T_circ) — the AIBECS bundled loaders
# already strip units inside ext/AIBECSJLD2Ext.jl, so T is a plain
# SparseMatrixCSC{Float64}.
# ---------------------------------------------------------------------------

const TIERS = (
    small  = [("OCIM0", OCIM0.load),
              ("OCCA",  OCCA.load)],
    medium = [("OCIM1", OCIM1.load),
              ("OCIM2", OCIM2.load)],
    large  = [("OCIM2_48L", OCIM2_48L.load)],
)
