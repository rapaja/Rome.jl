# Romeo LTI systems
#
# System type hierarchy.

"""
    Base type representing all SISO LTI systems.
"""
abstract type SisoLtiSystem end

"""
    Zero-system
    
A system that blocks all input signals: ideal no-pass system.
"""
struct ZeroSys <: SisoLtiSystem end

"""
    Unit-system

A system that passes all input signals: ideal all-pass system.
"""
struct UnitSys <: SisoLtiSystem end

"""
    Differintegral.

Differentiator of order `α` (when α>0), integrator of order `-α` (when α<0), or
identity operator (when α=0).

Note that when the order `α` is actually zero, a special object of type `UnitSys`
is returned.
"""
struct Diff{T <: Number} <: SisoLtiSystem
    α::T

    Diff{T}(α) where {T <: Number} = (α ≠ 0) ? new(α) : UnitSys()
end

Diff(α::T) where {T <: Number} = Diff{T}(α)

"""
    A SISO LTI system obtained by series connection of a scaler (pure gain) and another SISO LTI system.
"""
struct ScaledSystem{T <: Number} <: SisoLtiSystem
    k::T
    inner::SisoLtiSystem

    function ScaledSystem{T}(k, inner) where {T <: Number}
        if k == 1
            inner
        elseif k == 0
            ZeroSys()
        else
            new(k, inner)
        end # if
    end # function
end

ScaledSystem(k::T, inner::SisoLtiSystem) where {T <: Number} = ScaledSystem{T}(k, inner)
ScaledSystem(k::T, inner::ScaledSystem) where {T <: Number} = ScaledSystem{T}(k*inner.k, inner.inner)

"""
    A SISO LTI system obtained by series connection of two other SISO LTI systems.
"""
struct SeriesSystem <: SisoLtiSystem
    first::SisoLtiSystem
    second::SisoLtiSystem

    SeriesSystem(first::SisoLtiSystem, second::SisoLtiSystem) = new(first, second)
    SeriesSystem(first::SisoLtiSystem, second::UnitSys) = first
    SeriesSystem(first::UnitSys, second::SisoLtiSystem) = second
    SeriesSystem(first::SisoLtiSystem, second::ZeroSys) = ZeroSys()
    SeriesSystem(first::ZeroSys, second::SisoLtiSystem) = ZeroSys()
    SeriesSystem(first::SisoLtiSystem, second::ScaledSystem{<:Number}) = ScaledSystem(second.k, first*second.inner)
    SeriesSystem(first::ScaledSystem{<:Number}, second::SisoLtiSystem) = ScaledSystem(first.k, first.inner*second)
    SeriesSystem(first::ScaledSystem{<:Number}, second::ScaledSystem{<:Number}) = ScaledSystem(first.k*second.k, first.inner*second.inner)
end


"""
    A SISO LTI system obtained by parallel connection of two other SISO LTI systems.
"""
struct ParallelSystem <: SisoLtiSystem
    first::SisoLtiSystem
    second::SisoLtiSystem

    function ParallelSystem(first::SisoLtiSystem, second::SisoLtiSystem)
        if first == second
            2 * first
        else
            new(first, second)
        end
    end

    ParallelSystem(first::SisoLtiSystem, second::ZeroSys) = first
    ParallelSystem(first::ZeroSys, second::SisoLtiSystem) = second
    ParallelSystem(first::SisoLtiSystem, second::UnitSys) = first + 1
    ParallelSystem(first::UnitSys, second::SisoLtiSystem) = second + 1

    function ParallelSystem(first::SisoLtiSystem, second::ScaledSystem{<:Number})
        if second.inner == first
            (second.k+1) * first
        else
            new(first, second)
        end
    end

    ParallelSystem(first::ScaledSystem{<:Number}, second::SisoLtiSystem) = ParallelSystem(second, first)

    function ParallelSystem(first::ScaledSystem{<:Number}, second::ScaledSystem{<:Number})
        if first.inner == second.inner
            (first.k + second.k) * first.inner
        else
            new(first, second)
        end
    end
end

"""
    A SISO LTI system obtained by dividing two SISO LTI systems in the complex domain.
"""
struct RationalSystem <: SisoLtiSystem
    num::SisoLtiSystem
    den::SisoLtiSystem
end

"""
    A SISO LTI system obtained by taking the square root of another system in the complex domain.
"""
struct SquareRootSystem <: SisoLtiSystem
    inner::SisoLtiSystem
end

"""
    A SISO LTI system obtained by taking a power of another system in the complex domain.
"""
struct PowerSystem{T <: Number} <: SisoLtiSystem
    α::T
    inner::SisoLtiSystem
end


Base.:+(sys::SisoLtiSystem) = sys

Base.:-(sys::SisoLtiSystem) = ScaledSystem(-1, sys)
Base.:-(sys::ScaledSystem) = ScaledSystem(-sys.k, sys.inner)

Base.:*(first::SisoLtiSystem, second::SisoLtiSystem) = SeriesSystem(first, second)
Base.:*(num::Number, sys::SisoLtiSystem) = ScaledSystem(num, sys)
Base.:*(sys::SisoLtiSystem, num::Number) = num * sys

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

Base.sqrt(sys::SisoLtiSystem) = SquareRootSystem(sys)
Base.sqrt(sys::PowerSystem) = PowerSystem(sys.α/2, sys.inner)

Base.:^(sys::Diff{T}, n::Number) where T <: Number = Diff{typeof(one(T) * n)}(sys.α * n)

function Base.:^(sys::SquareRootSystem, n::Integer)
    if n == 0
        return Diff(0)
    elseif  n == 1
        return sys
    end
    m = abs(n)
    G = iseven(m) ? (sys.inner)^(m/2) : sqrt(sys.inner) * (sys.inner)^((m-1)/2)
    n>0 ? G : 1/G
end

function Base.:^(sys::SquareRootSystem, n::Number)
    if n == 0
        return Diff(0)
    elseif  n == 1
        return sys
    end
    if trunc(Int, n) == n
        return sys.^trunc(Int, n)
    else
        return PowerSystem(n/2, sys.inner)
    end
end

function Base.:^(sys::PowerSystem, n::Number)
    if n == 0
        return Diff(0)
    elseif n == 1
        return sys
    end
    α = sys.α * n;
    if trunc(Int, α) == α
        return sys.inner^trunc(Int, α)
    elseif trunc(Int, 2*α) == 2*α
        return SquareRootSystem(sys.inner)^trunc(Int, 2*α)
    else
        return PowerSystem(α, sys.inner)
    end
end

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

Base.:^(sys::SisoLtiSystem, n::Number) = PowerSystem(n, sys)

export SisoLtiSystem, Diff, ScaledSystem, ParallelSystem, SeriesSystem, RationalSystem, SquareRootSystem, PowerSystem