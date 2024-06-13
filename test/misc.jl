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

    @test A.n_couplings(model) == 2
    cbs = A.coupling_bounds(model)
    @test cbs in [([-4.321, -5], [1.234, -5]), ([-5, -4.321], [-5, 1.234])]
    @test merge(A.coupling_weights.(Ref(model), A.couplings(model))...) == Dict(
        "EX_o2_e" => 1,
        "EX_co2_e" => -1,
        "EX_glc__D_e" => 1,
        "EX_fru_e" => 1,
        "EX_pyr_e" => 1,
    )
end

@testset "Ugly cases of coupling" begin
    m = M.MATFBCModel(
        "model",
        Dict(
            "S" => [
                1.0 1.0
                -1.0 1.0
            ],
            "b" => [3.0, 1.0],
            "c" => [-0.25, 1.0],
            "xl" => -ones(2),
            "xu" => 2.0 * ones(2),
            "rxns" => ["r$x" for x = 1:2],
            "mets" => ["m$x" for x = 1:2],
        ),
    )

    @test size(A.coupling(m)) == (0, 2)

    mc = deepcopy(m)

    mc.mat["C"] = [0.5 0.5]
    @test A.coupling_bounds(mc) == ([-Inf], [Inf])

    m1 = deepcopy(mc)
    m1.mat["A"] = vcat(m1.mat["S"], m1.mat["C"])
    m1.mat["b"] = vcat(m1.mat["b"], [5])
    m1.mat["csense"] = ["L"]
    @test A.coupling_bounds(m1) == ([-Inf], [5.0])

    m2 = deepcopy(mc)
    m2.mat["d"] = [5]
    m2.mat["dsense"] = ["G"]
    @test A.coupling_bounds(m2) == ([5.0], [Inf])
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
