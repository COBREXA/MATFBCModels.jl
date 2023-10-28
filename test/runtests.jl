using Test
import AbstractFBCModels as A
import MATFBCModels as M

isdir("downloaded") || mkdir("downloaded")

iml1515_path = A.download_data_file(
    "http://bigg.ucsd.edu/static/models/iML1515.mat",
    joinpath("downloaded", "iML1515.mat"),
    "223db0b1ed69a7cc782f2a3093c1f48911a3c8bbd5bdf4bcdb13185cab5fdaa0",
)

iml1515_path = joinpath("test", "downloaded", "iML1515.mat")

@testset "MATFBCModels" begin

end
