@testset "Miscellaneous tests" begin

    model = A.load(M.MATFBCModel, joinpath(@__DIR__, "test-models", "iML1515.mat"))

    @test all(
        x in
        ["description" "c" "rev" "mets" "grRules" "subSystems" "b" "metFormulas" "rxnGeneMat" "S" "metNames" "lb" "metCharge" "ub" "rxnNames" "rxns" "genes"]
        for x in keys(model.mat)
    )

    @test all(x in ["mat" "MAT"] for x in A.filename_extensions(M.MATFBCModel))
end


@testset "Test coupling" begin
    model = convert(
        A.CanonicalModel.Model,
        A.load(M.MATFBCModel, joinpath(@__DIR__, "test-models", "e_coli_core.mat")),
    )
    model.couplings["breathing"] = A.CanonicalModel.Coupling(
        lower_bound = -4.321,
        upper_bound = 1.234,
        reaction_weights = Dict("EX_o2_e" => 1, "EX_co2_e" => -1),
    )
    model.couplings["eating"] = A.CanonicalModel.Coupling(
        lower_bound = -5,
        upper_bound = -5,
        reaction_weights = Dict("EX_glc__D_e" => 1, "EX_fru_e" => 1, "EX_pyr_e" => 1),
    )
    path = joinpath(@__DIR__, "test-models", "e_coli_core_coupled.mat")
    A.save(convert(M.MATFBCModel, model), path)
    A.run_fbcmodel_file_tests(
        M.MATFBCModel,
        path;
        name = "E. coli coupled",
        test_save = true,
    )
    model = convert(A.CanonicalModel.Model, A.load(M.MATFBCModel, path))
    # TODO @test if the couplings were carried over properly (except for the
    # IDs, these are discarded)
    # TODO add entries to `model.mat` manually to test if the other 2 ways of
    # expressing the couplings work
end

@testset "Corner cases" begin
    import MATFBCModels:
        parse_charge, eval_gene_association, flatten_gene_association, sortunique

    @test parse_charge(1) == 1
    @test parse_charge(2.0) == 2
    @test parse_charge("3") == 3
    @test parse_charge("4.0") == 4
    @test parse_charge(nothing) == nothing
    @test_throws ArgumentError parse_charge("totally positive charge")
    @test_throws DomainError parse_charge(["very charged"])

    @test_throws DomainError eval_gene_association(:(xor(gene("a"), gene("b"))), _ -> false)
    @test_throws DomainError flatten_gene_association(:(xor(gene("a"), gene("b"))))
    @test sortunique([3, 2, 2, 1]) == [1, 2, 3]
end
