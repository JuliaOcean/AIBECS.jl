
#=============================================
    Generate Dummy user-supplied model setup
=============================================#

using SparseArrays, LinearAlgebra, SuiteSparse

# model sizes
const nb = 3 # number of boxes
const m = 4 # number of parameters
const nt = 5 # number of tracers
const n = nb * nt

# assign local sources minus sinks to each tracer
for t in 1:nt
    fun = Symbol(string("G", t))
    vars = [Symbol(string("x", t2)) for t2 in 1:nt]
    push!(vars, :p)
    @eval $fun($(vars...)) = $t * $(vars[mod1(t+1, nt)]) .- $(:p)[mod1($t, m)] * $(vars[mod1(t-1, nt)])
end
Gs = tuple([@eval($(Symbol(string("G", t)))) for t in 1:nt]...)

# assign a circulation to each tracer
T0 = sprand(nb, nb, 0.5)
for t in 1:nt
    fun = Symbol(string("T", t))
    @eval $fun(p) = $(:p)[mod1($t, m)] * $T0 + 10 * $(:p)[mod1($t, m)] * I
end
Ts = tuple([@eval($(Symbol(string("T", t)))) for t in 1:nt]...)

# parameters
p = exp.(randn(m))
# random state vector
x = rand(n)

# observations (+ volume)
μx = [exp(randn()) * x[j:j+nb-1] for j in 1:nb:nb*nt]
σ²x = [exp(randn()) * x[j:j+nb-1] / 10 for j in 1:nb:nb*nt]
μp = [exp(randn()) * p[j] for j in 1:m]
σ²p = [exp(randn()) * p[j] / 10 for j in 1:m]
v = 1e3 * exp.(randn(nb))
ωs = tuple(exp.(randn(nt))...)
ωp = exp(rand())

#=============================================
    Generate Dummy user-supplied model setup
=============================================#

F, ∇ₓF = build_F_and_∇ₓF(Ts, Gs, nt, nb)
f, ∇ₓf = build_f_and_∇ₓf(ωs, μx, σ²x, v, ωp, μp, σ²p)
