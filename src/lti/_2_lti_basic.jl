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


export dkernel_identity, dkernel_difference