"""
$(TYPEDEF)

Wrapper around the models loaded in dictionaries from the MATLAB representation.

# Fields
$(TYPEDFIELDS)
"""
struct MATFBCModel <: A.AbstractFBCModel
    id::String
    mat::Dict{String,Any}
end

A.n_metabolites(m::MATFBCModel) = length(m.mat[guesskeys(:metabolites, m)])::Int
A.n_reactions(m::MATFBCModel) = length(m.mat[guesskeys(:reactions, m)])::Int
A.n_genes(m::MATFBCModel) = length(m.mat[guesskeys(:genes, m)])::Int

A.reactions(m::MATFBCModel) =
    String[x for x in vec(m.mat[guesskeys(:reactions, m)])]::Vector{String}
A.metabolites(m::MATFBCModel) =
    String[x for x in vec(m.mat[guesskeys(:metabolites, m)])]::Vector{String}
A.genes(m::MATFBCModel) =
    String[x for x in vec(m.mat[guesskeys(:genes, m)])]::Vector{String}

A.stoichiometry(m::MATFBCModel) =
    sparse(m.mat[guesskeys(:stoichiometry, m)])::SparseMatrixCSC{Float64,Int64}

A.bounds(m::MATFBCModel) = (
    Float64[x for x in vec(m.mat[guesskeys(:lbs, m)])],
    Float64[x for x in vec(m.mat[guesskeys(:ubs, m)])],
)::Tuple{Vector{Float64},Vector{Float64}}

A.balance(m::MATFBCModel) =
    sparsevec(m.mat[guesskeys(:balance, m)])::SparseVector{Float64,Int64}

A.objective(m::MATFBCModel) =
    sparsevec(m.mat[guesskeys(:objective, m)])::SparseVector{Float64,Int64}

A.reaction_gene_products_available(model::MATFBCModel, rid::String, available::Function) =
    A.reaction_gene_products_available_from_dnf(model, rid, available)

function A.reaction_gene_association_dnf(m::MATFBCModel, rid::String)
    if any(haskey(m.mat, x) for x in constants.keynames.grrs)
        grr = m.mat[guesskeys(:grrs, m)][findfirst(==(rid), A.reactions(m))]
        typeof(grr) == String ? parse_grr(grr) : nothing
    else
        nothing
    end
end

function A.metabolite_formula(m::MATFBCModel, mid::String)
    idx = findfirst(==(mid), A.metabolites(m))
    parse_formula(m.mat[guesskeys(:metabolite_formulas, m)][idx])
end

function A.metabolite_charge(m::MATFBCModel, mid::String)
    idx = findfirst(==(mid), A.metabolites(m))
    parse_charge(m.mat[guesskeys(:metabolite_charges, m)][idx])
end

A.metabolite_compartment(m::MATFBCModel, mid::String) = nothing

function A.reaction_stoichiometry(m::MATFBCModel, rid::String)
    ridx = first(indexin([rid], m.mat[guesskeys(:reactions, m)]))[1] # get the index out of the cartesian index
    met_inds = findall(m.mat[guesskeys(:stoichiometry, m)][:, ridx] .!= 0.0)
    Dict{String,Float64}(
        m.mat[guesskeys(:metabolites, m)][met_ind] =>
            m.mat[guesskeys(:stoichiometry, m)][met_ind, ridx] for met_ind in met_inds
    )
end

function A.reaction_name(m::MATFBCModel, rid::String)
    idx = findfirst(==(rid), A.reactions(m))
    isnothing(idx) && return nothing
    string(m.mat[guesskeys(:reaction_names, m)][idx])::String
end

function A.metabolite_name(m::MATFBCModel, mid::String)
    idx = findfirst(==(mid), A.metabolites(m))
    isnothing(idx) && return nothing
    string(m.mat[guesskeys(:metabolite_names, m)][idx])::String
end

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

function Base.convert(::Type{MATFBCModel}, m::A.AbstractFBCModel)
    if typeof(m) == MATFBCModel
        return m
    end

    lb, ub = A.bounds(m)
    return MATFBCModel(
        "model", # default name
        Dict(
            "S" => A.stoichiometry(m),
            "rxns" => A.reactions(m),
            "mets" => A.metabolites(m),
            "lb" => Vector(lb),
            "ub" => Vector(ub),
            "b" => Vector(A.balance(m)),
            "c" => Vector(A.objective(m)),
            "genes" => A.genes(m),
            "grRules" => [
                unparse_grr(A.reaction_gene_association_dnf(m, rid)) for
                rid in A.reactions(m)
            ],
            "metFormulas" =>
                [unparse_formula(A.metabolite_formula(m, mid)) for mid in A.metabolites(m)],
            "metCharges" => [A.metabolite_charge(m, mid) for mid in A.metabolites(m)],
        ),
    )
end
