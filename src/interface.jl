"""
$(TYPEDEF)

Wrapper around the models loaded in dictionaries from the MATLAB representation.

# Fields
$(TYPEDFIELDS)
"""
struct MATFBCModel <: A.AbstractFBCModel
    mat::Dict{String,Any}
end

A.n_metabolites(m::MATFBCModel) = size(m.mat["S"], 1)
A.n_reactions(m::MATFBCModel) = size(m.mat["S"], 2)
A.n_genes(m::MATFBCModel) = length(m.mat["genes"])

function A.reactions(m::MATFBCModel)
    if haskey(m.mat, "rxns")
        reshape(m.mat["rxns"], n_reactions(m))
    else
        "rxn" .* string.(1:n_reactions(m))
    end
end

function A.metabolites(m::MATFBCModel)
    nm = n_metabolites(m)
    if haskey(m.mat, "mets")
        reshape(m.mat["mets"], length(m.mat["mets"]))[begin:nm]
    else
        "met" .* string.(1:n_metabolites(m))
    end
end

A.stoichiometry(m::MATFBCModel) = sparse(m.mat["S"])

A.bounds(m::MATFBCModel) = (
    reshape(get(m.mat, "lb", fill(-Inf, n_reactions(m), 1)), n_reactions(m)),
    reshape(get(m.mat, "ub", fill(Inf, n_reactions(m), 1)), n_reactions(m)),
)

function A.balance(m::MATFBCModel)
    b = get(m.mat, "b", spzeros(n_metabolites(m), 1))
    sparse(reshape(b, n_metabolites(m)))
end

A.objective(m::MATFBCModel) =
    sparse(reshape(get(m.mat, "c", zeros(n_reactions(m), 1)), n_reactions(m)))

function A.genes(m::MATFBCModel)
    x = get(m.mat, "genes", [])
    reshape(x, length(x))
end

function A.reaction_gene_associations(m::MATFBCModel, rid::String)
    if haskey(m.mat, "grRules")
        grr = m.mat["grRules"][findfirst(==(rid), reactions(m))]
        typeof(grr) == String ? A.parse_grr(grr) : nothing
    else
        nothing
    end
end

A.metabolite_formula(m::MATFBCModel, mid::String) = A.maybemap(
    x -> A.parse_formula(x[findfirst(==(mid), metabolites(m))]),
    A.gets(m.mat, nothing, Internal.constants.keynames.metformulas),
)

function A.metabolite_charge(m::MATFBCModel, mid::String)
    met_charge = A.maybemap(
        x -> x[findfirst(==(mid), metabolites(m))],
        A.gets(m.mat, nothing, Internal.constants.keynames.metcharges),
    )
    A.maybemap(Int, isnan(met_charge) ? nothing : met_charge)
end

function A.metabolite_compartment(m::MATFBCModel, mid::String)
    res = A.maybemap(
        x -> x[findfirst(==(mid), metabolites(m))],
        A.gets(m.mat, nothing, Internal.constants.keynames.metcompartments),
    )
    # if the metabolite is an integer or a (very integerish) float, it is an
    # index to a table of metabolite names (such as in the yeast GEM)
    typeof(res) <: Real || return res
    return A.maybemap(
        table -> table[Int(res)],
        A.gets(m.mat, nothing, Internal.constants.keynames.metcomptables),
    )
end

function A.reaction_stoichiometry(m::MATFBCModel, rid::String)
    ridx = first(indexin([rid], m.mat["rxns"]))[1] # get the index out of the cartesian index
    met_inds = findall(m.mat["S"][:, ridx] .!= 0.0)
    Dict(m.mat["mets"][met_ind] => m.mat["S"][met_ind, ridx] for met_ind in met_inds)
end

A.reaction_name(m::MATFBCModel, rid::String) = A.maybemap(
    x -> x[findfirst(==(rid), reactions(m))],
    A.gets(m.mat, nothing, Internal.constants.keynames.rxnnames),
)

A.metabolite_name(m::MATFBCModel, mid::String) = A.maybemap(
    x -> x[findfirst(==(mid), metabolites(m))],
    A.gets(m.mat, nothing, Internal.constants.keynames.metnames),
)

A.gene_name(m::MATFBCModel, gid::String) = nothing

# NOTE: There's no useful standard on how and where to store notes and
# annotations in MATLAB models. We therefore leave it very open for the users,
# who can easily support any annotation scheme using a custom wrapper.
# Even the (simple) assumptions about grRules, formulas and charges that we use
# here are very likely completely incompatible with >50% of the MATLAB models
# out there.
A.gene_annotations(model::MATFBCModel, gid::String) = A.Annotations()
A.gene_notes(model::MATFBCModel, gid::String) = A.Notes()
A.reaction_annotations(model::MATFBCModel, rid::String) = A.Annotations()
A.reaction_notes(model::MATFBCModel, rid::String) = A.Notes()
A.metabolite_annotations(model::MATFBCModel, mid::String) = A.Annotations()
A.metabolite_notes(model::MATFBCModel, mid::String) = A.Notes()

# function Base.convert(::Type{MATFBCModel}, m::AbstractMetabolicModel)
#     if typeof(m) == MATFBCModel
#         return m
#     end

#     lb, ub = variable_bounds(m)
#     cl, cu = coupling_bounds(m)
#     nr = n_reactions(m)
#     nm = n_metabolites(m)
#     return MATFBCModel(
#         Dict(
#             "S" => stoichiometry(m),
#             "rxns" => reactions(m),
#             "mets" => metabolites(m),
#             "lb" => Vector(lb),
#             "ub" => Vector(ub),
#             "b" => Vector(metabolite_bounds(m)),
#             "c" => Vector(objective(m)),
#             "C" => coupling(m),
#             "cl" => Vector(cl),
#             "cu" => Vector(cu),
#             "genes" => genes(m),
#             "grRules" =>
#                 default.(
#                     "",
#                     maybemap.(
#                         x -> unparse_grr(String, x),
#                         reaction_gene_associations.(Ref(m), reactions(m)),
#                     ),
#                 ),
#             "metFormulas" =>
#                 default.(
#                     "",
#                     maybemap.(
#                         unparse_formula,
#                         metabolite_formula.(Ref(m), metabolites(m)),
#                     ),
#                 ),
#             "metCharges" => default.(0, metabolite_charge.(Ref(m), metabolites(m))),
#             "metCompartments" =>
#                 default.("", metabolite_compartment.(Ref(m), metabolites(m))),
#         ),
#     )
# end
