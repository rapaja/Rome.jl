"""
    convolve(g, u)

Convolution of signals `g` and `u`.
"""
function convolve(g::Vector{T}, u::Vector{T}) where T <: Number
    @assert length(g) == length(u) "Romeo.LTI.Operators.convolve: Input signals must be of the same length."
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
Compute signal `g` such that `g ⋆ u = y`.
"""
function deconvolve(y::Vector{T}, u::Vector{T}) where T <: Number
    @assert length(y) == length(u) "Romeo.LTI.Operators.deconvolve: Input signals must be of the same length."
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

"""
    convroot(g)

Convolutional root of `h`.
Compute signal `g` such that `h = g ⋆ g`.

For signals having negative initial values, the convolutive root is multivalued
(similar to the algebraic root). The convolutive root will return the solution
obtained using the primary branch of the algebraic root.
"""
function convroot(h::Vector{T}) where T <: Number
    if length(h) == 0
        g = Vector{typeof(one(T) / one(T))}(undef, size(h))
    else
        @assert h[1] != 0 "Romeo.LTI.Operators.convroot: The initial value of the inpur signal cannot be zero."
        if T <: Real
            if h[1] > 0
                g = Vector{typeof(one(T) / one(T))}(undef, size(h))
            else
                g = Vector{typeof(one(Complex{T}) / one(Complex{T}))}(undef, size(h))
            end
        else
            g = Vector{typeof(one(T) / one(T))}(undef, size(h))
        end
        if length(g) > 0
            g[1] = sqrt(h[1])
            if length(g) >= 2
                @inbounds @simd for k = 2:length(h)
                    val = h[k]
                    @simd for i = 2:(k - 1)
                        val -= g[k - i + 1] * g[i]
                    end
                    g[k] = val / 2 / g[1]
                end
            end
        end
    end
    return g
end

"""
    convid([T=typeof(1.0)], N)

Convolutional identity operator.

Kernel length is `N`. Element type is determined by `T`.
"""
function convid(T::Type{<:Number}, N::Integer)
    δ = zeros(T, N)
    δ[1] = one(T)
    return δ
end # function

convid(N::Integer) = convid(typeof(1.0), N)

"""
    convinv(g)

Convolutional inverse of the signal `g`.
"""
function convinv(g::Vector{T}) where T <: Number
    id = convid(T, length(g))
    return deconvolve(id, g)
end # end
    
export convolve, deconvolve, convroot, convid, convinv