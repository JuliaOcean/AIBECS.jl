module AIBECSRecipesParametersExt

using AIBECS
using RecipesBase
using Distributions
using Unitful: ustrip, upreferred
using UnPack

@userplot PlotParameter
@recipe function f(plt::PlotParameter)
    d, v, initv, s, u = extract_dvisu(plt.args...)
    xs = default_range(d)
    ys = [pdf(d, x) for x in xs]
    xu = xs * upreferred(u) .|> u .|> ustrip
    vu = v * upreferred(u) |> u |> ustrip
    @series begin
        label --> "prior"
        xu, ys
    end
    @series begin
        label --> "initial value"
        initv * [1,1], [0, pdf(d, v)]
    end
    @series begin
        xguide --> "$s ($(string(u)))"
        yguide --> "PDF"
        label --> "value"
        markershape --> :o
        [vu], [pdf(d, v)]
    end
end

@userplot PlotParameters
@recipe function f(plt::PlotParameters)
    p, = plt.args
    layout --> (length(p),1)
    for (i,s) in enumerate(AIBECS.flattenable_symbols(p))
        d, v, initvu, s, u = extract_dvisu(p,s)
        xs = default_range(d)
        ys = [pdf(d, x) for x in xs]
        xu = xs * upreferred(u) .|> u .|> ustrip
        vu = v * upreferred(u) |> u |> ustrip
        initv = ustrip(upreferred(initvu * AIBECS.units(p,s)))
        @series begin
            label --> "prior"
            subplot := i
            xu, ys
        end
        @series begin
            label --> "initial value"
            subplot := i
            initvu * [1,1], [0, pdf(d, initv)]
        end
        @series begin
            xguide --> "$s ($(string(u)))"
            yguide --> "PDF"
            label --> "value"
            markershape --> :o
            subplot := i
            [vu], [pdf(d, v)]
        end
    end
end

extract_dvisu(d, v, initv, s, u) = d, v, initv, s, u
extract_dvisu(p::AIBECS.AbstractParameters, s) = AIBECS.prior(p,s), UnPack.unpack(p, Val(s)), AIBECS.initial_value(p,s), s, AIBECS.units(p,s)

function default_range(d::Distribution, n = 4)
    μ, σ = mean(d), std(d)
    xmin = max(μ - n*σ, minimum(support(d)))
    xmax = min(μ + n*σ, maximum(support(d)))
    range(xmin, xmax, length=1001)
end

end # module
