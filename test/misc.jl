@testset "Miscellaneous tests" begin

    model = A.load(M.MATFBCModel, joinpath(@__DIR__, "test-models", "iML1515.mat"))

    @test all(
        in.(
            keys(model.mat),
            Ref(
                [
                    "description"
                    "c"
                    "rev"
                    "mets"
                    "grRules"
                    "subSystems"
                    "b"
                    "metFormulas"
                    "rxnGeneMat"
                    "S"
                    "metNames"
                    "lb"
                    "metCharge"
                    "ub"
                    "rxnNames"
                    "rxns"
                    "genes"
                ],
            ),
        ),
    )

    @test all(in.(A.filename_extensions(M.MATFBCModel), Ref([
        "mat"
        "MAT"
    ])))
end
