
"""
    module key_names

This module collects information on how people typically name stuff in the MAT
file structure.
"""
module key_names
reactions = ["reactions", "rxns", "RXNS", "REACTIONS", "Reactions", "Rxns"]

metabolites = ["metabolites", "mets", "METS", "METABOLITES", "Metabolites", "Mets"]

genes = ["genes", "GENES", "Genes"]

lbs = ["lbs", "lb", "lowerbounds", "lower_bounds"]

ubs = ["ubs", "ub", "upperbounds", "upper_bounds"]

stoichiometry = ["S"]

balance = ["b"]

objective = ["c"]

grrs = ["gene_reaction_rules", "grRules", "rules"]

ids = ["id", "description"]

metabolite_formulas = ["metFormula", "metFormulas"]

metabolite_charges = ["metCharge", "metCharges"]

reaction_names = ["rxnNames"]

metabolite_names = ["metNames"]
end
