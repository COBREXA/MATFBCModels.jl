
"""
$(TYPEDEF)

A wrapper around a matlab model loaded in the MATLAB representation as defined
by `MAT.jl` package.

# Fields
$(TYPEDFIELDS)
"""
struct MATFBCModel <: A.AbstractFBCModel
    id::String
    mat::Dict{String,Any}
end

const Maybe = A.Maybe
