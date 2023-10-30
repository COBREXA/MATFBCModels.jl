@testset "Miscellaneous tests" begin

    model = A.load(M.MATFBCModel, joinpath(@__DIR__, "test-models", "iML1515.mat"))

    @test all(
        x in
        ["description" "c" "rev" "mets" "grRules" "subSystems" "b" "metFormulas" "rxnGeneMat" "S" "metNames" "lb" "metCharge" "ub" "rxnNames" "rxns" "genes"]
        for x in keys(model.mat)
    )

    @test all(x in ["mat" "MAT"] for x in A.filename_extensions(M.MATFBCModel))
end
@testset "Corner cases" begin
    import MATFBCModels: parse_charge

    @test parse_charge(1) == 1
    @test parse_charge(2.0) == 2
    @test parse_charge("3") == 3
    @test parse_charge("4.0") == 4
    @test parse_charge(nothing) == nothing
    @test_throws ArgumentError parse_charge("totally positive charge")
    @test_throws DomainError parse_charge(["very charged"])
end
