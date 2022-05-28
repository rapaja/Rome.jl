# Romeo LTI systems
#
# System type hierarchy.

# TODO: Remove `SquareRootSystem` from hierarchy!
#
#   This system type is special only regarding kernel operation.
#   It can be handled as special case there.
#   It only increases complexity in the rest of the code.

# TODO: Parallel and series type systems should be able to have multiple inner systems instead of just two.

# TODO: Consider introducing `PolynomialSystem` as a special case of `ParallelSystem`.

"""
    Base type representing all SISO LTI systems.
"""
abstract type SisoLtiSystem end

Base.:+(sys::SisoLtiSystem) = sys
Base.:-(sys::SisoLtiSystem) = -1 * sys


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
struct Diff{T<:Number} <: SisoLtiSystem
    α::T

    Diff{T}(α) where {T<:Number} = (α ≠ 0) ? new(α) : UnitSys()
end

Diff(α::T) where {T<:Number} = Diff{T}(α)

Base.:(==)(first::Diff{<:Number}, second::Diff{<:Number}) = first.α == second.α

"""
    A SISO LTI system obtained by series connection of a scaler (pure gain) and another SISO LTI system.
"""
struct ScaledSystem{T<:Number} <: SisoLtiSystem
    k::T
    inner::SisoLtiSystem

    function ScaledSystem{T}(k, inner) where {T<:Number}
        if k == 1
            inner
        elseif k == 0 || inner == ZeroSys()
            ZeroSys()
        else
            new(k, inner)
        end # if
    end # function
end

ScaledSystem(k::T, inner::SisoLtiSystem) where {T<:Number} = ScaledSystem{T}(k, inner)
ScaledSystem(k::T, inner::ScaledSystem{<:Number}) where {T<:Number} = ScaledSystem{T}(k * inner.k, inner.inner)

Base.:(==)(first::ScaledSystem{<:Number}, second::ScaledSystem{<:Number}) =
    first.inner == second.inner && first.k == second.k

Base.:*(num::Number, sys::SisoLtiSystem) = ScaledSystem(num, sys)
Base.:*(sys::SisoLtiSystem, num::Number) = num * sys

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
    SeriesSystem(first::SisoLtiSystem, second::ScaledSystem{<:Number}) = ScaledSystem(second.k, first * second.inner)
    SeriesSystem(first::ScaledSystem{<:Number}, second::SisoLtiSystem) = ScaledSystem(first.k, first.inner * second)
    SeriesSystem(first::ScaledSystem{<:Number}, second::ScaledSystem{<:Number}) = ScaledSystem(first.k * second.k, first.inner * second.inner)
    SeriesSystem(first::Diff{<:Number}, second::Diff{<:Number}) = Diff(first.α + second.α)
end

Base.:*(first::SisoLtiSystem, second::SisoLtiSystem) = SeriesSystem(first, second)

"""
    A SISO LTI system obtained by parallel connection of two other SISO LTI systems.
"""
struct ParallelSystem <: SisoLtiSystem
    first::SisoLtiSystem
    second::SisoLtiSystem

    ParallelSystem(first::SisoLtiSystem, second::SisoLtiSystem) = (first == second) ? 2 * first : new(first, second)
    ParallelSystem(first::SisoLtiSystem, second::ZeroSys) = first
    ParallelSystem(first::ZeroSys, second::SisoLtiSystem) = second
    ParallelSystem(first::SisoLtiSystem, second::ScaledSystem{<:Number}) = (second.inner == first) ? (second.k + 1) * first : new(first, second)
    ParallelSystem(first::ScaledSystem{<:Number}, second::SisoLtiSystem) = ParallelSystem(second, first)
    ParallelSystem(first::ScaledSystem{<:Number}, second::ScaledSystem{<:Number}) = (first.inner == second.inner) ? (first.k + second.k) * first.inner : new(first, second)
end

Base.:+(first::SisoLtiSystem, second::SisoLtiSystem) = ParallelSystem(first, second)
Base.:+(sys::SisoLtiSystem, num::Number) = sys + num * Diff(0)
Base.:+(num::Number, sys::SisoLtiSystem) = sys + num

Base.:-(first::SisoLtiSystem, second::SisoLtiSystem) = first + (-second)
Base.:-(sys::SisoLtiSystem, num::Number) = sys + (-num)
Base.:-(num::Number, sys::SisoLtiSystem) = num + (-sys)

"""
    A SISO LTI system obtained by taking the square root of another system in the complex domain.
"""
struct SquareRootSystem <: SisoLtiSystem
    inner::SisoLtiSystem

    SquareRootSystem(inner::SisoLtiSystem) = new(inner)
    SquareRootSystem(inner::ZeroSys) = ZeroSys()
    SquareRootSystem(inner::UnitSys) = UnitSys()
    SquareRootSystem(inner::Diff{<:Number}) = Diff(inner.α / 2)
    SquareRootSystem(inner::ScaledSystem{<:Number}) = ScaledSystem(sqrt(inner.k), sqrt(inner.inner))
end

Base.sqrt(sys::SisoLtiSystem) = SquareRootSystem(sys)

"""
    A SISO LTI system obtained by taking a power of another system in the complex domain.
"""
struct PowerSystem{T<:Number} <: SisoLtiSystem
    α::T
    inner::SisoLtiSystem

    function PowerSystem{T}(α::T, inner::SisoLtiSystem) where {T<:Number}
        if α == 0
            UnitSys()
        elseif α == 1
            inner
        elseif α == 0.5
            SquareRootSystem(inner)
        else
            new(α, inner)
        end # if
    end # function

    PowerSystem{T}(α::T, inner::ZeroSys) where {T<:Number} = (α ≠ 0) ? ZeroSys() : UnitSys()
    PowerSystem{T}(α::T, inner::UnitSys) where {T<:Number} = UnitSys()
    PowerSystem{T}(α::T, inner::Diff{<:Number}) where {T<:Number} = Diff(inner.α * α)
    PowerSystem{T}(α::T, inner::ScaledSystem{<:Number}) where {T<:Number} = ScaledSystem(inner.k^α, inner.inner^α)
    PowerSystem{T}(α::T, inner::SquareRootSystem) where {T<:Number} = PowerSystem(α / 2, inner.inner)
end

PowerSystem(α::T, inner::SisoLtiSystem) where {T<:Number} = PowerSystem{T}(α, inner)

SquareRootSystem(inner::SquareRootSystem) = PowerSystem(1 // 4, inner)
SquareRootSystem(inner::PowerSystem{<:Number}) = PowerSystem(inner.α / 2, inner.inner)

Base.:^(sys::SisoLtiSystem, n::Number) = PowerSystem(n, sys)

"""
    A SISO LTI system obtained by dividing two SISO LTI systems in the complex domain.
"""
struct RationalSystem <: SisoLtiSystem
    num::SisoLtiSystem
    den::SisoLtiSystem

    function RationalSystem(num::SisoLtiSystem, den::SisoLtiSystem)
        if num == den
            UnitSys()
        else
            new(num, den)
        end
    end

    RationalSystem(num::SisoLtiSystem, den::UnitSys) = num
    RationalSystem(num::UnitSys, den::SisoLtiSystem) = PowerSystem(-1, den)
    RationalSystem(num::SisoLtiSystem, den::ZeroSys) = @error "Cannot create rational systems with zero denominator."
    RationalSystem(num::ZeroSys, den::SisoLtiSystem) = ZeroSys()
    RationalSystem(num::SisoLtiSystem, den::ScaledSystem{<:Number}) = ScaledSystem(1 / den.k, num / den.inner)
    RationalSystem(num::ScaledSystem{<:Number}, den::SisoLtiSystem) = ScaledSystem(num.k, num.inner / den)
    RationalSystem(num::ScaledSystem{<:Number}, den::ScaledSystem{<:Number}) = ScaledSystem(num.k / den.k, num.inner / den.inner)
    RationalSystem(num::Diff{<:Number}, den::Diff{<:Number}) = Diff(num.α - den.α)
    RationalSystem(num::SisoLtiSystem, den::PowerSystem{<:Number}) = (den.inner == num) ? PowerSystem(1 - den.α, den.inner) : new(num, den)
    RationalSystem(num::PowerSystem{<:Number}, den::SisoLtiSystem) = (num.inner == den) ? PowerSystem(num.α - 1, num.inner) : new(num, den)
    RationalSystem(num::PowerSystem{<:Number}, den::PowerSystem{<:Number}) = (num.inner == den.inner) ? PowerSystem(num.α - den.α, num.inner) : new(num, den)
end

Base.:/(num::SisoLtiSystem, den::SisoLtiSystem) = RationalSystem(num, den)
Base.:/(sys::SisoLtiSystem, den::Number) = ScaledSystem(1 / den, sys)
Base.:/(num::Number, sys::SisoLtiSystem) = RationalSystem(num * Diff(0), sys)

Base.:\(den::SisoLtiSystem, num::SisoLtiSystem) = RationalSystem(num, den)
Base.:\(sys::SisoLtiSystem, num::Number) = RationalSystem(num * Diff(0), sys)
Base.:\(den::Number, sys::SisoLtiSystem) = ScaledSystem(1 / den, sys)

export SisoLtiSystem, Diff, ScaledSystem, ParallelSystem, SeriesSystem, RationalSystem, SquareRootSystem, PowerSystem