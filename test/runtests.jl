using Test
import RequiredInterfaces as R
import AbstractFBCModels as A
import MATFBCModels as M

isdir("downloaded") || mkdir("downloaded")

iml1515_path = A.download_data_file(
    "http://bigg.ucsd.edu/static/models/iML1515.mat",
    joinpath("downloaded", "iML1515.mat"),
    "223db0b1ed69a7cc782f2a3093c1f48911a3c8bbd5bdf4bcdb13185cab5fdaa0",
)

@testset "MATFBCModels" begin
    @testset "Interface implemented" R.check_implementations(A.AbstractFBCModel)
    # IO
    include("io.jl")
    # accessors
    include("accessors.jl")
end
