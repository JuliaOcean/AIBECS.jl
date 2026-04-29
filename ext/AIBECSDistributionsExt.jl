module AIBECSDistributionsExt

using AIBECS
using Distributions: ContinuousUnivariateDistribution, logpdf
using Bijectors: transformed

AIBECS.mismatch(d::ContinuousUnivariateDistribution, λ) = -logpdf(transformed(d), λ)

end # module
