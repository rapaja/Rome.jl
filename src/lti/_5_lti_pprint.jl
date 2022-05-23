# Romeo LTI pretty printing.
#
# Utility functions for human-readable representations of systems.

"""
    pprint(sys [, var="s"])

Generates pretty human-readable representation of Laplace transform of the system `sys`,
with the Laplace variable denoted by `var`. 
"""
pprint(sys::SisoLtiSystem) = pprint(sys, 's')

function pprint(sys::Diff{<: Number}, var::Char)
    if sys.α == 0
        return "1"
    elseif sys.α == 1
        return "$var"
    else
        return "$var^$(sys.α)"
    end
end

pprint(sys::ScaledSystem, var::Char) = "$(sys.k) * ($(pprint(sys.inner, var)))"
pprint(sys::SeriesSystem, var::Char) = "($(pprint(sys.first, var))) * ($(pprint(sys.second, var)))"
pprint(sys::ParallelSystem, var::Char) = "$(pprint(sys.first, var)) + $(pprint(sys.second, var))"
pprint(sys::RationalSystem, var::Char) = "($(pprint(sys.first, var))) / ($(pprint(sys.second, var)))"
pprint(sys::SquareRootSystem, var::Char) = "sqrt($(pprint(sys.inner, var)))"
pprint(sys::PowerSystem, var::Char) = "($(pprint(sys.inner, var)))^$(sys.α)"

export pprint