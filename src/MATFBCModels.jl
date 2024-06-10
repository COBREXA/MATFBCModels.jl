module MATFBCModels

using DocStringExtensions

import AbstractFBCModels as A

using MAT
using SparseArrays

include("types.jl")
include("constants.jl")
include("interface.jl")
include("io.jl")
include("utils.jl")
include("grr_utils.jl")

export MATFBCModel

end # module MatFBCModels
