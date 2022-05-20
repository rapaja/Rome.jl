import SpecialFunctions

"""
    Base type representing all SISO LTI systems.
"""
abstract type SisoLtiSystem end

"""
    kernel([T=typeof(1.0), ] sys, N, Δt)

Given a system `sys` compute its kernel vector of length `N`, uniformly
discretized with sample time `Δt`.

When simulating the response `y` of any system to a given signal `u`, the kernel is 
defined here as a sequence of samples `k` such that the response is equal to the 
discrete convolution `k ⋆ u`.
"""
kernel(sys::SisoLtiSystem, N::Integer, Δt::Real) = kernel(typeof(1.0), sys, N, Δt)

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
    simulate(sys, input, Δt)

Simulate response of `sys` to input signal `input`, uniformly sampled with `Δt`.
"""
function simulate(sys::SisoLtiSystem, input::Array{T}, Δt::Real) where T <: Number
    convolve(kernel(T, sys, length(input), Δt), input)
end

"""
Differintegral.

Differentiator of order `α` (when α>0), integrator of order `-α` (when α<0), or
identity operator (when α=0).
"""
struct Diff{T <: Number} <: SisoLtiSystem
    α::T
end

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
    A SISO LTI system obtained by series connection of a scaler (pure gain) and another SISO LTI system.
"""
struct ScaledSystem{T <: Number} <: SisoLtiSystem
    k::T
    inner::SisoLtiSystem
end

"""
    A SISO LTI system obtained by parallel connection of two other SISO LTI systems.
"""
struct ParallelSystem <: SisoLtiSystem
    first::SisoLtiSystem
    second::SisoLtiSystem
end

"""
    A SISO LTI system obtained by series connection of two other SISO LTI systems.
"""
struct SeriesSystem <: SisoLtiSystem
    first::SisoLtiSystem
    second::SisoLtiSystem
end

"""
    A SISO LTI system obtained by dividing two SISO LTI systems.
"""
struct RationalSystem <: SisoLtiSystem
    num::SisoLtiSystem
    den::SisoLtiSystem
end

function kernel(T::Type{<:Number}, sys::Diff{<: Number}, N::Integer, Δt::Real)
    if sys.α == 0
        res = zeros(T, N)
        res[1] = one(T)
        return res
    elseif sys.α == 1
        res = zeros(T, N)
        res[1] = one(T)
        res[2] = -one(T)
        return res / Δt
    elseif real(sys.α) > 0
        ddt = kernel(T, Diff(one(typeof(sys.α))), N, Δt)
        temp = kernel(T, Diff(sys.α - 1), N, Δt)
        return convolve(ddt, temp)
    elseif real(sys.α) < 0
        ik = t -> t^(-sys.α) * SpecialFunctions.gamma(-sys.α + 1)
        t = 0:Δt:(N * Δt)
        ikt = ik.(t)
        return ikt[2:end] - ikt[1:end - 1]
    else
        @error "Not implemented!"
    end
end


kernel(T::Type{<:Number}, sys::ScaledSystem, N::Integer, Δt::Real) = sys.k * kernel(T, sys.inner, N, Δt)
kernel(T::Type{<:Number}, sys::ParallelSystem, N::Integer, Δt::Real) = kernel(T, sys.first, N, Δt) .+ kernel(T, sys.second, N, Δt)
kernel(T::Type{<:Number}, sys::SeriesSystem, N::Integer, Δt::Real) = convolve(kernel(T, sys.first, N, Δt), kernel(T, sys.second, N, Δt))
kernel(T::Type{<:Number}, sys::RationalSystem, N::Integer, Δt::Real) = deconvolve(kernel(T, sys.num, N, Δt), kernel(T, sys.den, N, Δt))

Base.:+(sys::SisoLtiSystem) = sys
Base.:-(sys::SisoLtiSystem) = ScaledSystem(-1, sys)

Base.:*(first::SisoLtiSystem, second::SisoLtiSystem) = SeriesSystem(first, second)
Base.:*(sys::SisoLtiSystem, num::Number) = ScaledSystem(num, sys)
Base.:*(num::Number, sys::SisoLtiSystem) = ScaledSystem(num, sys)

Base.:+(first::SisoLtiSystem, second::SisoLtiSystem) = ParallelSystem(first, second)
Base.:+(sys::SisoLtiSystem, num::Number) = sys + num * Diff(0)
Base.:+(num::Number, sys::SisoLtiSystem) = sys + num * Diff(0)

Base.:-(first::SisoLtiSystem, second::SisoLtiSystem) = first + (-second)
Base.:-(sys::SisoLtiSystem, num::Number) = sys + (-num)
Base.:-(num::Number, sys::SisoLtiSystem) = num + (- sys)

Base.:/(num::SisoLtiSystem, den::SisoLtiSystem) = RationalSystem(num, den)
Base.:/(sys::SisoLtiSystem, den::Number) = ScaledSystem(1 / den, sys)
Base.:/(num::Number, sys::SisoLtiSystem) = RationalSystem(num * Diff(0), sys)

Base.:\(den::SisoLtiSystem, num::SisoLtiSystem) = RationalSystem(num, den)
Base.:\(sys::SisoLtiSystem, num::Number) = RationalSystem(num * Diff(0), sys)
Base.:\(den::Number, sys::SisoLtiSystem) = ScaledSystem(1 / den, sys)

function Base.:^(sys::SisoLtiSystem, n::Integer)
    if n == 0
        return Diff(0)
    elseif n == 1
        return sys
    elseif n > 1
        return sys * sys^(n - 1)
    else # negative `n`
        return (1 / sys)^(-n)
    end
end

Base.:^(sys::Diff{T}, n::Integer) where T <: Number = Diff{typeof(one(T) * n)}(sys.α * n)
Base.:^(sys::Diff{T}, n::Number) where T <: Number = Diff{typeof(one(T) * n)}(sys.α * n)


export SisoLtiSystem, Diff, ScaledSystem, ParallelSystem, SeriesSystem, RationalSystem
export kernel, kernel_to_signal, impulse_resp, step_resp, simulate, invlaplace