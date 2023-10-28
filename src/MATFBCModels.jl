module MATFBCModels

using DocStringExtensions
import AbstractFBCModels as A
using MAT, SparseArrays

include("constants.jl")
include("interface.jl")
include("io.jl")

export MATFBCModel

end # module MatFBCModels
