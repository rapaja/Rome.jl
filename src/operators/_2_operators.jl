"""
    dkernel_identity([T=typeof(1.0)], N)

Diskrete kernel of the **identity operator**.

Kernel length is `N`. Element type is determined by `T`.
"""
function dkernel_identity(T::Type{<:Number}, N::Integer)
    dkernel = zeros(T, N)
    dkernel[1] = one(T)
    return dkernel
end # function

dkernel_identity(N::Integer) = dkernel_identity(typeof(1.0), N)

"""
    dkernel_difference([T=typeof(1.0)], N)

Discrete kernel (of length `N`) of the **first order backward difference operator**.

Kernel length is `N`. Element type is determined by `T`.
"""
function dkernel_difference(T::Type{<:Number}, N::Integer)
    dkernel = zeros(T, N)
    dkernel[1] = one(T)
    dkernel[2] = -one(T)
    return dkernel
end # function

dkernel_difference(N::Integer) = dkernel_difference(typeof(1.0), N)

"""
    dkernel_inverse(g::Vector{T} where T <: Number)

Inverse kernel for the given one.
"""
function dkernel_inverse(g::Vector{T} where T <: Number)
    id = dkernel_identity(eltype(g), length(g))
    return deconvolve(id, g)
end # end

"""
    response(g, u)

Given the discrete kernel `g` and the input signal `u`,
compute covolution `g ⋆ u`.
"""
    response = convolve

    export response, dkernel_identity, dkernel_difference