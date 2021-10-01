"""
    convolve(g, u)

Convolution of signals `g` and `u`.
"""
function convolve(g::Vector{T}, u::Vector{T})::Vector{T} where T <: Number
    @assert length(g) == length(u) "Romeo.Operators.convolve: Input signals must be of the same length."
    y = similar(u)
    if length(y) > 0
        @inbounds @simd for t = 1:length(u)
            val = zero(y[t])
            @simd for τ = 1:t
                val += g[t - τ + 1] * u[τ]
            end # for
            y[t] = val
        end # for
    end # if
    return y
end # function

"""
    deconvolve(y, u)

Deconvolution of signals `y` and `u`.
"""
function deconvolve(y::Vector{T}, u::Vector{T}) where T <: Number
    @assert length(y) == length(u) "Romeo.Operators.deconvolve: Input signals must be of the same length."
    g = Vector{typeof(one(T) / one(T))}(undef, size(u))
    if length(u) > 0
        m = one(T) / u[1]
        g[1] = y[1] * m
        @inbounds @simd for t = 2:length(u)
            val = y[t]
            @simd for τ = 2:t
                val -= g[t - τ + 1] * u[τ]
            end # for
            g[t] = val * m
        end # for
    end # if
    return g
end # function

export convolve