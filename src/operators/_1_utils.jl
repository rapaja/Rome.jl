"""
    convolve(g, u)

Convolution of signals `g` and `u`.
"""
function convolve(g::Vector{T}, u::Vector{T})::Vector{T} where T <: Number
    y = similar(u)
    for t = 1:length(u)
        val = zero(y[t])
        for τ = 1:t
            val += g[t - τ + 1] * u[τ]
        end # for
        y[t] = val
    end # for
    return y
end

"""
    deconvolve(y, u)

Deconvolution of signals `y` and `u`.
"""
function deconvolve(y::Vector{T}, u::Vector{T})::Vector{T} where T <: Number
    g = similar(u)
    g[1] = y[1] / u[1]
    for t = 2:length(u)
        val = y[t]
        for τ = 2:t
            val -= g[t - τ + 1] * u[τ]
        end # for
        g[t] = val / u[1]
    end # for
    return g
end # function

export convolve