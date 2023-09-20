@testset "Accessors" begin
    model = M.load(M.MATFBCModel, iml1515_path)

    @test "SHK3Dr" in M.reactions(model)
    @test M.n_reactions(model) == 2712  
    @test "co2_e" in M.metabolites(model)
    @test M.n_metabolites(model) == 1877 
    @test "b0008" in M.genes(model) 
    @test M.n_genes(model) == 1516
    @test all(size(M.stoichiometry(model)) .== (1877, 2712))
    @test all(length.(M.bounds(model)) .== 2712)
    @test all(M.balance(model) .== 0) 
    @test M.objective(model)[2669] == 1
    @test all(in.(M.reaction_gene_associations(model, "FBA"), Ref([ ["b2925"],["b2097"]])))
    @test M.metabolite_formula(model, "atp_c")["C"] == 10
    @test M.metabolite_charge(model, "atp_c") == -4 
    @test isnothing(M.metabolite_compartment(model, "atp_c"))
    @test isempty(M.gene_annotations(model, "b0008"))
    @test isempty(M.gene_notes(model, "b0008"))
    @test isempty(M.reaction_annotations(model, "FBA"))
    @test isempty(M.reaction_notes(model, "FBA"))
    @test isempty(M.metabolite_annotations(model, "atp_c"))
    @test isempty(M.metabolite_notes(model, "atp_c"))
    @test isnothing(M.gene_name(model, "b0008"))
    @test M.reaction_stoichiometry(model, "FBA")["dhap_c"] == 1.0 
    @test M.reaction_name(model, "FBA") == "Fructose-bisphosphate aldolase"
    @test M.metabolite_name(model, "atp_c") == "ATP C10H12N5O13P3"
end
