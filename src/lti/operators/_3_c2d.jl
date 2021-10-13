"""
    dkernel_from_intkernelfcn(intkernelfcn, t)

Compute discrete convolution kernel given the integrated kernel function.

The discrete kernel is computed using the _backward Euler_ (i.e. _left rectangular_)
discretization. Therefore, the first element of the discrete kernel is always zero.
"""
function dkernel_from_intkernelfcn(intkernelfcn::Function, t::Array{T})::Array{T}  where T <: Real
    dkernel = similar(t)
    dkernel[1] = zero(T)
    dkernel[2:end] = intkernelfcn(t[2:end]) - intkernelfcn(t[1:end - 1])
    return dkernel
end

# function kernel_c2d(g, t)
#     dt = similar(t)
#     dt[1:end - 1] = diff(t)
#     dt[end] = dt[end - 1]
#     return g(t) .* dt
# end # function

export dkernel_from_intkernelfcn