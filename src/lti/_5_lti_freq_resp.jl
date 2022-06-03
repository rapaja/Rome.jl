# Romeo LTI frequency responses.
#
# Implementation of `freqresp` function of various systems.

freq_resp(sys::Diff{<:Number}, ω::AbstractVector{<:Number}) = (1im * ω) .^ sys.α

export freq_resp