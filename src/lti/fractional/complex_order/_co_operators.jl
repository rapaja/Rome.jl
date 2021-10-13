using SpecialFunctions

function kernel_fcn_co_integral(t, α, β)
    gmm = gamma(α + β * 1im)
    a, b = real(gmm), imag(gmm)
    lnt = log.(t)
    return t.^α .* (a * cos.(lnt) .- b * sin.(lnt))
end # function

export kernel_fcn_co_integral