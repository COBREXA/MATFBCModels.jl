using Test
import AbstractFBCModels as A
import MATFBCModels as M
import SparseArrays

@testset "MATFBCModels" begin

    A.run_fbcmodel_type_tests(M.MATFBCModel)

    modeldir = joinpath(@__DIR__, "test-models")
    mkpath(modeldir)

    for (name, url, hash, ts) in [
        (
            "e_coli_core",
            "http://bigg.ucsd.edu/static/models/e_coli_core.mat",
            "478e6fa047ede1b248975d7565208ac9363a44dd64aad1900b63127394f4175b",
            true,
        ),
        (
            "iJO1366",
            "http://bigg.ucsd.edu/static/models/iJO1366.mat",
            "b5cfe21b6369a00e45d600b783f89521f5cc953e25ee52c5f1d0a3f83743be30",
            true,
        ),
        (
            "iML1515",
            "http://bigg.ucsd.edu/static/models/iML1515.mat",
            "223db0b1ed69a7cc782f2a3093c1f48911a3c8bbd5bdf4bcdb13185cab5fdaa0",
            true,
        ),
    ]
        path = joinpath(modeldir, "$name.mat")
        A.download_data_file(url, path, hash)
        A.run_fbcmodel_file_tests(M.MATFBCModel, path; name, test_save = ts)
    end

    include("test_iML1515.jl")
    include("misc.jl")
end
