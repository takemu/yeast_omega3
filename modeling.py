import os
import pandas as pd
import cobra
from cobra import Reaction, Metabolite
import reframed
import ast
import re

evo_base_constr = {'R_r_1580': (0, 1000), 
                   'R_r_0488': (0, 1000), 
                   'R_r_0472': (0, 1000), 
                   'R_r_0438': (0, 1000), 
                   'R_r_1021': (-1000, 1000), 
                   'R_r_0226': (0, 1000), 
                   'R_r_0662': (0, 1000), 
                   'R_r_0502': (-1000, 1000), 
                   'R_r_0503': (-1000, 1000), 
                   'R_r_0659': (-1000, 1000), 
                   'R_r_0732': (-1000, 1000), 
                   'R_r_0733': (-1000, 1000)}
reframed_to_gpr_rxns = {}  # {reframed_rxn: [gpr_rxn1, gpr_rxn2, ...]}
new_reframed_rxn_dict = {}  # {new_reframed_rxn: [r_name, reversible]}


def add_reaction(rxn, r_name, lb, ub, metabolites_to_add, gene_reaction_rule=None):
    reaction = Reaction(rxn, r_name, lower_bound=lb, upper_bound=ub)
    reaction.add_metabolites(metabolites_to_add)
    if gene_reaction_rule:
        reaction.gene_reaction_rule = gene_reaction_rule
    return reaction


def add_omega3_pathway(cobra_model):
    s_o001 = Metabolite(id='s_o001', name='DHA', formula='C22H32O2', charge=0, compartment='c')
    s_o002 = Metabolite(id='s_o002', name='DHA', formula='C22H32O2', charge=0, compartment='e')
    cobra_model.add_metabolites([s_o001, s_o002])    
    
    s_1260 = cobra_model.metabolites.get_by_id('s_1260') # oleate
    s_1101 = cobra_model.metabolites.get_by_id('s_1101') # malonyl-CoA_c
    s_1212 = cobra_model.metabolites.get_by_id('s_1212') # NADPH_c
    s_0794 = cobra_model.metabolites.get_by_id('s_0794') # H+_c
    s_1275 = cobra_model.metabolites.get_by_id('s_1275') # O2_c
    s_0529 = cobra_model.metabolites.get_by_id('s_0529') # CoA_c
    s_1203 = cobra_model.metabolites.get_by_id('s_1203') # NADH_c
    s_0803 = cobra_model.metabolites.get_by_id('s_0803') # H2O_c

    gene1 = cobra.Gene('DHAS')  # artificial gene for DHA synthase
    gene1.name = 'DHA synthase'
    cobra_model.genes.append(gene1)
    r_o001 = add_reaction('r_o001', 'DHA synthesis', 0, 1000, 
                          {s_1260: -1.0, s_1101: -4.0, s_1212: -5.0, s_0794: -5.0, s_1275: -5.0, s_o001: 1.0, s_0529: 4.0, s_1203: 5.0, s_0803: 4.0}, 'DHAS')
    r_o002 = add_reaction('r_o002', 'DHA transport', -1000, 1000, {s_o001: -1.0, s_o002: 1.0})
    r_o003 = add_reaction('r_o003', 'DHA exchange', 0, 1000, {s_o002: -1.0})
    cobra_model.add_reactions([r_o001, r_o002, r_o003])
    return cobra_model

def create_gpr_model(model_file):
    reframed_model_file = f'{os.path.dirname(model_file)}/reframed-GEM-omega3.xml'

    if not os.path.exists(reframed_model_file):
        cobra_model = cobra.io.read_sbml_model(model_file)
        cobra_model = add_omega3_pathway(cobra_model)
        cobra.io.write_sbml_model(cobra_model, reframed_model_file)

    reframed_model = reframed.load_cbmodel(reframed_model_file, flavor='fbc2')
    gpr_model = reframed.core.transformation.gpr_transform(reframed_model, inplace=False, add_proteome=True, gene_prefix='G_', usage_prefix='u_', pseudo_genes=None)

    gpr_rxn_df = pd.DataFrame([[r[0], r[1].name] for r in gpr_model.reactions.items()], columns=['id', 'r_name']).set_index('id')
    for reframed_rxn, _ in reframed_model.reactions.items():
        reframed_to_gpr_rxns[reframed_rxn] = gpr_rxn_df[gpr_rxn_df.index.str.startswith(reframed_rxn)].index.tolist()
    '''Create a Dataframe indexed by isozyme-distinguished reframed reactions'''
    for gpr_rxn, row in gpr_rxn_df.iterrows():
        new_reframed_rxn = gpr_rxn
        direction = ''
        m = re.match(r'(R\_r\_[^_]+)(.*)$', gpr_rxn)
        if m:
            new_reframed_rxn = m.group(1)
            direction = m.group(2)
            m2 = re.match(r'(R\_r\_[^_]+)(.*)(\_iso\d+)$', gpr_rxn)
            if m2:
                new_reframed_rxn += m2.group(3)
                direction = m2.group(2)
        reversible = False
        if direction in ['_f', '_b']:
            reversible = True
        if new_reframed_rxn in new_reframed_rxn_dict:
            new_reframed_rxn_dict[new_reframed_rxn][1] = reversible
        else:
            new_reframed_rxn_dict[new_reframed_rxn] = [row['r_name'], reversible]

    return gpr_model, reframed_model
