abstract type LTISystem end

kernel(sys::LTISystem, N::Integer) = kernel(Float64, sys, N)

kernel_to_signal(k::Array{<:Number}, h::Real) = k / h

function simulate(sys::LTISystem, input::Array{T}, h::Real) where T <: Number
    kernel_to_signal(convolve(kernel(T, sys, length(input)), input), h)
end

struct Diff{T <: Number} <: LTISystem
    α::T
end

impulse_resp(T::Type{<:Number}, sys::LTISystem, N::Integer, h::Real) = kernel_to_signal(kernel(T, sys, N), h)
impulse_resp(sys::LTISystem, N::Integer, h::Real) = impulse_resp(Float64, sys::LTISystem, N::Integer, h::Real)

step_resp(T::Type{<:Number}, sys::LTISystem, N::Integer, h::Real) = impulse_resp(T, sys / (Diff(1) / h), N, h)
step_resp(sys::LTISystem, N::Integer, h::Real) = step_resp(Float64, sys, N, h)

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

function kernel(T::Type{<:Number}, sys::Diff{<: Number}, N::Integer)
    if sys.α == 0
        res = zeros(T, N)
        res[1] = one(T)
        return res
    elseif sys.α == 1
        res = zeros(T, N)
        res[1] = one(T)
        res[2] = -one(T)
        return res
    else
        @error "Not implemented!"
    end
end

kernel(T::Type{<:Number}, sys::ScaledSystem, N::Integer) = sys.k * kernel(T, sys.inner, N)
kernel(T::Type{<:Number}, sys::ParallelSystem, N::Integer) = kernel(T, sys.first, N) .+ kernel(T, sys.second, N)
kernel(T::Type{<:Number}, sys::SeriesSystem, N::Integer) = convolve(kernel(T, sys.first, N), kernel(T, sys.second, N))
kernel(T::Type{<:Number}, sys::RationalSystem, N::Integer) = deconvolve(kernel(T, sys.num, N), kernel(T, sys.den, N))

unit = Diff(0)

Base.:+(sys::LTISystem) = sys
Base.:-(sys::LTISystem) = ScaledSystem(-1, sys)

Base.:*(first::LTISystem, second::LTISystem) = SeriesSystem(first, second)
Base.:*(sys::LTISystem, num::Number) = ScaledSystem(num, sys)
Base.:*(num::Number, sys::LTISystem) = ScaledSystem(num, sys)

Base.:+(first::LTISystem, second::LTISystem) = ParallelSystem(first, second)
Base.:+(sys::LTISystem, num::Number) = sys + num * unit
Base.:+(num::Number, sys::LTISystem) = sys + num * unit

Base.:-(first::LTISystem, second::LTISystem) = first + (-second)
Base.:-(sys::LTISystem, num::Number) = sys + (-num)
Base.:-(num::Number, sys::LTISystem) = num + (- sys)

Base.:/(num::LTISystem, den::LTISystem) = RationalSystem(num, den)
Base.:/(sys::LTISystem, den::Number) = ScaledSystem(1 / den, sys)
Base.:/(num::Number, sys::LTISystem) = RationalSystem(num * unit, sys)

Base.:\(den::LTISystem, num::LTISystem) = RationalSystem(num, den)
Base.:\(sys::LTISystem, num::Number) = RationalSystem(num * unit, sys)
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


export LTISystem, Diff, ScaledSystem, ParallelSystem, SeriesSystem, RationalSystem
export kernel, kernel_to_signal, impulse_resp, step_resp