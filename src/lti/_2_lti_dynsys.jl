import SpecialFunctions

abstract type LTISystem end

kernel(sys::LTISystem, N::Integer, Δt::Real) = kernel(Float64, sys, N, Δt)

kernel_to_signal(k::Array{<:Number}, Δt::Real) = k / Δt

function simulate(sys::LTISystem, input::Array{T}, Δt::Real) where T <: Number
    convolve(kernel(T, sys, length(input), Δt), input)
end

struct Diff{T <: Number} <: LTISystem
    α::T
end

impulse_resp(T::Type{<:Number}, sys::LTISystem, N::Integer, Δt::Real) = kernel_to_signal(kernel(T, sys, N, Δt), Δt)
impulse_resp(sys::LTISystem, N::Integer, Δt::Real) = impulse_resp(Float64, sys::LTISystem, N::Integer, Δt::Real)

invlaplace(T::Type{<:Number}, sys::LTISystem, N::Integer, Δt::Real) = impulse_resp(T, sys, N, Δt)
invlaplace(sys::LTISystem, N::Integer, Δt::Real) = impulse_resp(sys, N, Δt)

step_resp(T::Type{<:Number}, sys::LTISystem, N::Integer, Δt::Real) = impulse_resp(T, sys / Diff(1), N, Δt)
step_resp(sys::LTISystem, N::Integer, Δt::Real) = step_resp(Float64, sys, N, Δt)

struct ScaledSystem{T <: Number} <: LTISystem
    k::T
    inner::LTISystem
end

struct ParallelSystem <: LTISystem
    first::LTISystem
    second::LTISystem
end

struct SeriesSystem <: LTISystem
    first::LTISystem
    second::LTISystem
end

struct RationalSystem <: LTISystem
    num::LTISystem
    den::LTISystem
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

Base.:+(sys::LTISystem) = sys
Base.:-(sys::LTISystem) = ScaledSystem(-1, sys)

Base.:*(first::LTISystem, second::LTISystem) = SeriesSystem(first, second)
Base.:*(sys::LTISystem, num::Number) = ScaledSystem(num, sys)
Base.:*(num::Number, sys::LTISystem) = ScaledSystem(num, sys)

Base.:+(first::LTISystem, second::LTISystem) = ParallelSystem(first, second)
Base.:+(sys::LTISystem, num::Number) = sys + num * Diff(0)
Base.:+(num::Number, sys::LTISystem) = sys + num * Diff(0)

Base.:-(first::LTISystem, second::LTISystem) = first + (-second)
Base.:-(sys::LTISystem, num::Number) = sys + (-num)
Base.:-(num::Number, sys::LTISystem) = num + (- sys)

Base.:/(num::LTISystem, den::LTISystem) = RationalSystem(num, den)
Base.:/(sys::LTISystem, den::Number) = ScaledSystem(1 / den, sys)
Base.:/(num::Number, sys::LTISystem) = RationalSystem(num * Diff(0), sys)

Base.:\(den::LTISystem, num::LTISystem) = RationalSystem(num, den)
Base.:\(sys::LTISystem, num::Number) = RationalSystem(num * Diff(0), sys)
Base.:\(den::Number, sys::LTISystem) = ScaledSystem(1 / den, sys)

function Base.:^(sys::LTISystem, n::Integer)
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


export LTISystem, Diff, ScaledSystem, ParallelSystem, SeriesSystem, RationalSystem
export kernel, kernel_to_signal, impulse_resp, step_resp, simulate, invlaplace