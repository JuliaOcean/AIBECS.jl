# transportMatrixTests.jl


# Age test, should be ≈ 1000 years
function computeIdealMeanAge(wet3d,T)
    l3d = copy(wet3d)
    l3d[:,:,2:end] .= 0
    iwet = find(wet3d)
  L = spdiagm(l3d[iwet])
    nwet = length(iwet)
  return (T + L) \ ones(nwet)
end
meanΓyrs = mean(computeIdealMeanAge(wet3d,T)) / (365 * spd)
@test 100 < meanΓyrs < 1000

# T 1 ≈ 0 test
Tcolmax = vec(maximum(abs.(T),1))
Trowmax = vec(maximum(abs.(T),2))
T1 = T * ones(nwet)
maximum(abs.(T1) ./ Tcolmax)
maximum(abs.(T1) ./ Trowmax)

# vᵀ T ≈ 0 test
vᵀT = T.' * v
maximum(abs.(vᵀT) ./ (v .* Tcolmax))
maximum(abs.(vᵀT) ./ (v .* Trowmax))

# S 1 ≈ 0 test
Scolmax = vec(maximum(abs.(S),1))
Srowmax = vec(maximum(abs.(S),2))
S1 = S * ones(nwet)
maximum(abs.(S1) ./ Scolmax)
maximum(abs.(S1) ./ Srowmax)
S13d = NaN * wet3d
S13d[iwet] .= S1
using Plots
gr()
plot(S13d[50,54,:],line = :scatter)

# vᵀ S ≈ 0 test
vᵀS = S.' * v
maximum(abs.(vᵀS) ./ (v .* Scolmax))
maximum(abs.(vᵀS) ./ (v .* Srowmax))

function isMassConserving(x,v)
  vᵀx = v.' * x
  return vᵀx
