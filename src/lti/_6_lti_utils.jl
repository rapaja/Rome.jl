# Romeo LTI utils
#
# Various utility functions.

"""
    kernel_to_signel(k, Δt)

Given a kernel vector `k` compute the corresponding response signal.

This function is an unfortunate necesity of the discretization process, and the consequence
of the way the continuous-time convolutions are implemented using sampled signals. Input
particular, if the response of a system to a given excitation `u` is computed as `g ⋆ u`,
with `⋆` denoting discrete convolution, then the continuous-time impulse response of the
system is not `g`, but `g/Δt`.
"""
kernel_to_signal(k::Array{<:Number}, Δt::Real) = k / Δt

"""
    impulse_resp([T=typeof(1.0)], sys, N, Δt)

Given a system `sys` compute its impulse response vector of length `N`, uniformly
discretized with sample time `Δt`.
"""
impulse_resp(T::Type{<:Number}, sys::SisoLtiSystem, N::Integer, Δt::Real) = kernel_to_signal(kernel(T, sys, N, Δt), Δt)
impulse_resp(sys::SisoLtiSystem, N::Integer, Δt::Real) = impulse_resp(typeof(1.0), sys::SisoLtiSystem, N::Integer, Δt::Real)

"""
    invlaplace([T=typeof(1.0)], sys, N, Δt)

Given a Laplace domain operator (a system) `sys` compute its inverse Laplace transform
as a vector of length `N`, uniformly discretized with sample time `Δt`.
"""
invlaplace(T::Type{<:Number}, sys::SisoLtiSystem, N::Integer, Δt::Real) = impulse_resp(T, sys, N, Δt)
invlaplace(sys::SisoLtiSystem, N::Integer, Δt::Real) = impulse_resp(sys, N, Δt)

"""
    step_resp([T=typeof(1.0)], sys, N, Δt)

Given a system `sys` compute its step response vector of length `N`, uniformly
discretized with sample time `Δt`.
"""
step_resp(T::Type{<:Number}, sys::SisoLtiSystem, N::Integer, Δt::Real) = impulse_resp(T, sys / Diff(1), N, Δt)
step_resp(sys::SisoLtiSystem, N::Integer, Δt::Real) = step_resp(Float64, sys, N, Δt)

"""
    simulate(sys, input, Δt)

Simulate response of `sys` to input signal `input`, uniformly sampled with `Δt`.
"""
function simulate(sys::SisoLtiSystem, input::Array{T}, Δt::Real) where T <: Number
    convolve(kernel(T, sys, length(input), Δt), input)
end


export kernel_to_signal, impulse_resp, step_resp, simulate, invlaplace