module MATFBCModels

using DocStringExtensions
import AbstractFBCModels as A
using MAT, SparseArrays, Reexport

module Internal
    using DocStringExtensions
    import ..A
    include("constants.jl")
end

import .Internal

include("model.jl")
include("io.jl")

export MATFBCModel

@reexport using AbstractFBCModels

end # module MatFBCModels
