
"""
A named tuple that contains the magic values that are used globally for
whatever purposes.
"""
const constants = (
    default_reaction_bound = 1e3,
    keynames = (
        rxns = ["reactions", "rxns", "RXNS", "REACTIONS", "Reactions", "Rxns"],
        mets = ["metabolites", "mets", "METS", "METABOLITES", "Metabolites", "Mets"],
        genes = ["genes", "GENES", "Genes"],
        lbs = ["lbs", "lb", "lowerbounds", "lower_bounds"],
        ubs = ["ubs", "ub", "upperbounds", "upper_bounds"],
        stochiometry = ["S"],
        balance = ["b"],
        objective = ["c"],
        grrs = ["gene_reaction_rules", "grRules", "rules"],
        ids = ["id", "description"],
        metformulas = ["metFormula", "metFormulas"],
        metcharges = ["metCharge", "metCharges"],
        metcompartments = ["metCompartment", "metCompartments", "metComp", "metComps"],
        metcomptables = ["comps", "compNames"],
        rxnnames = ["rxnNames"],
        metnames = ["metNames"],
    ),
)
