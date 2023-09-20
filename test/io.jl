@testset "IO" begin
    model = M.load(M.MATFBCModel, iml1515_path)

    @test all(
        in.(
            keys(model.mat),
            Ref([
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
            ]),
        ),
    )

    @test all(in.(M.filename_extensions(M.MATFBCModel), Ref([
        "mat"
        "MAT"
    ])))

    # TODO save

    # TODO convert
end
