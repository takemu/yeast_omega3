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

appl_base_constr = {'R_r_0487': (0, 0), # close glycerol utilization pathway as glucose is used as c-source
                   'R_r_0488': (0, 0), 
                   'R_r_0659': (0, 1000), # constrain isocitrate dehydrogenase in NADPH generating direction
                   'R_r_0502': (0, 1000), # glycine synthesis from serine when glucose is the carbon source
                   'R_r_0503': (0, 1000),
                   'R_r_0732': (0, 1000),
                   'R_r_0733': (0, 1000),
                   }

anaerobic_constr = {'R_r_0438': (0, 0),
                    'R_r_1021': (0, 0),
                    'R_r_0226': (0, 0),
                    'R_r_0662': (0, 0)}

reframed_to_gpr_rxns = {}  # {reframed_rxn: [gpr_rxn1, gpr_rxn2, ...]}
new_reframed_rxn_dict = {}  # {new_reframed_rxn: [r_name, reversible]}


def add_reaction(rxn, r_name, lb, ub, metabolites_to_add, gene_reaction_rule=None):
    reaction = Reaction(rxn, r_name, lower_bound=lb, upper_bound=ub)
    reaction.add_metabolites(metabolites_to_add)
    if gene_reaction_rule:
        reaction.gene_reaction_rule = gene_reaction_rule
    return reaction

def add_aromatic_pathway(cobra_model):
    # add_Ehrlich_yeast76.m
    s_a001 = Metabolite('s_a001', formula='C8H8O2', name='4-hydroxyphenylacetaldehyde', compartment='c')
    s_a002 = Metabolite('s_a002', formula='C8H8O3', name='4-hydroxyphenylacetate', compartment='c')
    s_a003 = Metabolite('s_a003', formula='C8H8O3', name='4-hydroxyphenylacetate', compartment='e')
    s_a004 = Metabolite('s_a004', formula='C9H8O4', name='4-hydroxyphenylpyruvate', compartment='c')
    s_a005 = Metabolite('s_a005', formula='C9H8O4', name='4-hydroxyphenylpyruvate', compartment='e')
    s_a006 = Metabolite('s_a006', formula='C8H10O2', name='tyrosol', compartment='c')
    s_a007 = Metabolite('s_a007', formula='C8H10O2', name='tyrosol', compartment='e')
    cobra_model.add_metabolites([s_a001, s_a002, s_a003, s_a004, s_a005, s_a006, s_a007])
    s_1051 = cobra_model.metabolites.get_by_id('s_1051')  # yeast76.mets{850}
    s_0180 = cobra_model.metabolites.get_by_id('s_0180')  # yeast76.mets{124}
    s_0991 = cobra_model.metabolites.get_by_id('s_0991')  # yeast76.mets{792}
    s_1399 = cobra_model.metabolites.get_by_id('s_1399')  # yeast76.mets{1088}
    s_0955 = cobra_model.metabolites.get_by_id('s_0955')  # yeast76.mets{756}
    s_0794 = cobra_model.metabolites.get_by_id('s_0794')  # yeast76.mets{605}
    s_0456 = cobra_model.metabolites.get_by_id('s_0456')  # yeast76.mets{343}
    s_0803 = cobra_model.metabolites.get_by_id('s_0803')  # yeast76.mets{614}
    s_1207 = cobra_model.metabolites.get_by_id('s_1207')  # yeast76.mets{961}
    s_1212 = cobra_model.metabolites.get_by_id('s_1212')  # yeast76.mets{965}
    s_1198 = cobra_model.metabolites.get_by_id('s_1198')  # yeast76.mets{952}
    s_1203 = cobra_model.metabolites.get_by_id('s_1203')  # yeast76.mets{957}

    r_a001 = add_reaction('r_a001', 'L-tyrosine transaminase (L-glu)', -1000, 1000, {s_1051: -1.0, s_0180: -1.0, s_a004: 1.0, s_0991: 1.0}, 'YGL202W')
    r_a002 = add_reaction('r_a002', 'L-tyrosine transaminase (L-ala)', -1000, 1000, {s_1051: -1.0, s_1399: -1.0, s_a004: 1.0, s_0955: 1.0}, 'YHR137W')
    r_a003 = add_reaction('r_a003', '4-hydroxyphenylpyruvate decarboxylase', 0, 1000, {s_a004: -1.0, s_0794: -1.0, s_0456: 1.0, s_a001: 1.0}, 'YGR087C or YLR044C or YLR134W')
    r_a004 = add_reaction('r_a004', 'aldehyde dehydrogenase (4-hydroxyphenylacetate)', -1000, 1000,
                            {s_a001: -1.0, s_0803: -1.0, s_1207: -1.0, s_a002: 1.0, s_0794: 1.0, s_1212: 1.0}, 'YPL061W')
    r_a005 = add_reaction('r_a005', 'aldehyde dehydrogenase (tyrosol)', -1000, 1000, {s_a001: -1.0, s_1198: 1.0, s_a006: 1.0, s_0794: -1.0, s_1203: -1.0},
                            'YBR145W or YDL168W or YOL086C')
    r_a006 = add_reaction('r_a006', 'tyrosol transport', -1000, 1000, {s_a007: 1.0, s_a006: -1.0})
    r_a007 = add_reaction('r_a007', 'tyrosol exchange', 0, 1000, {s_a007: -1.0})
    r_a008 = add_reaction('r_a008', '4-hydroxyphenylacetate transport', -1000, 1000, {s_a003: 1.0, s_a002: -1.0})
    r_a009 = add_reaction('r_a009', '4-hydroxyphenylacetate exchange', 0, 1000, {s_a003: -1.0})
    r_a010 = add_reaction('r_a010', '4-hydroxyphenylpyruvate transport', -1000, 1000, {s_a005: 1.0, s_a004: -1.0})
    r_a011 = add_reaction('r_a011', '4-hydroxyphenylpyruvate exchange', 0, 1000, {s_a005: -1.0})
    cobra_model.add_reactions([r_a001, r_a002, r_a003, r_a004, r_a005, r_a006, r_a007, r_a008, r_a009, r_a010, r_a011])
    return cobra_model

def add_indigoidine_pathway(cobra_model):
    s_i001 = Metabolite(id='s_i001', name='indigoidine', formula='C10H12N2O4', charge=0, compartment='c')
    s_i002 = Metabolite(id='s_i002', name='indigoidine', formula='C10H12N2O4', charge=0, compartment='e')
    cobra_model.add_metabolites([s_i001, s_i002])
    s_0999 = cobra_model.metabolites.get_by_id('s_0999')  # L-glutamine_c
    s_0434 = cobra_model.metabolites.get_by_id('s_0434')  # ATP_c
    s_0423 = cobra_model.metabolites.get_by_id('s_0423')  # AMP_c
    s_0633 = cobra_model.metabolites.get_by_id('s_0633')  # diphosphate_c
    s_0794 = cobra_model.metabolites.get_by_id('s_0794')  # H+_c

    gene1 = cobra.Gene('B8Y4G4')  # gene for indigoidine synthase
    gene1.name = 'bpsA'
    cobra_model.genes.append(gene1)

    r_i001 = add_reaction('r_i001', 'indigoidine synthase', 0, 1000, {s_0999: -2.0, s_0434: -1.0, s_i001: 1.0, s_0423: 1.0, s_0633: 1.0, s_0794: 2.0}, 'B8Y4G4')
    r_i002 = add_reaction('r_i002', 'indigoidine transport', -1000, 1000, {s_i001: -1.0, s_i002: 1.0})
    r_i003 = add_reaction('r_i003', 'indigoidine exchange', 0, 1000, {s_i002: -1.0})
    cobra_model.add_reactions([r_i001, r_i002, r_i003])
    return cobra_model

def add_bikaverin_pathway(cobra_model):
    s_b001 = Metabolite(id='s_b001', name='bikaverin', formula='C20H14O8', charge=0, compartment='c')
    s_b002 = Metabolite(id='s_b002', name='bikaverin', formula='C20H14O8', charge=0, compartment='e')
    cobra_model.add_metabolites([s_b001, s_b002])
    s_1101 = cobra_model.metabolites.get_by_id('s_1101')  # malonyl-CoA_c
    s_0373 = cobra_model.metabolites.get_by_id('s_0373')  # acetyl-CoA_c
    s_0529 = cobra_model.metabolites.get_by_id('s_0529')  # coenzyme A_c
    s_0456 = cobra_model.metabolites.get_by_id('s_0456')  # carbon dioxide_c

    gene1 = cobra.Gene('S0DZM7')  # gene for bikaverin synthase
    gene1.name = 'bik1'
    cobra_model.genes.append(gene1)

    r_b001 = add_reaction('r_b001', 'bikaverin synthase', 0, 1000, {s_1101: -7.0, s_0373: -1.0, s_b001: 1.0, s_0529: 8.0, s_0456: 7.0}, 'S0DZM7')
    r_b002 = add_reaction('r_b002', 'bikaverin transport', -1000, 1000, {s_b001: -1.0, s_b002: 1.0})
    r_b003 = add_reaction('r_b003', 'bikaverin exchange', 0, 1000, {s_b002: -1.0})
    cobra_model.add_reactions([r_b001, r_b002, r_b003])
    return cobra_model

def add_omega3_pathway(cobra_model):
    s_o001 = Metabolite(id='s_o001', name='DHA', formula='C22H32O2', charge=0, compartment='c')
    s_o002 = Metabolite(id='s_o002', name='DHA', formula='C22H32O2', charge=0, compartment='e')
    cobra_model.add_metabolites([s_o001, s_o002])    
    # s_1262 = cobra_model.metabolites.get_by_id('s_1262') # oleoyl-CoA_c
    s_1260 = cobra_model.metabolites.get_by_id('s_1260') # oleate
    s_1101 = cobra_model.metabolites.get_by_id('s_1101') # malonyl-CoA_c
    s_1212 = cobra_model.metabolites.get_by_id('s_1212') # NADPH_c
    s_0794 = cobra_model.metabolites.get_by_id('s_0794') # H+_c
    s_1275 = cobra_model.metabolites.get_by_id('s_1275') # O2_c
    s_0529 = cobra_model.metabolites.get_by_id('s_0529') # CoA_c
    s_1203 = cobra_model.metabolites.get_by_id('s_1203') # NADH_c
    s_0803 = cobra_model.metabolites.get_by_id('s_0803') # H2O_c

    # heterologous_genes = ['FAD2', 'FAD3', 'FADSD6', 'ELOVL5', 'FADSD5', 'ELOVL2', 'D4']
    # for gene in heterologous_genes:
    #     cobra_model.genes.append(cobra.Gene(gene))
    gene1 = cobra.Gene('DHAS')  # artificial gene for DHA synthase
    gene1.name = 'DHA synthase'
    cobra_model.genes.append(gene1)
    # r_o001 = add_reaction('r_o001', 'oleoyl-CoA hydrolysis', 0, 1000, {s_1262: -1.0, s_0803: -1.0, s_1260: 1.0, s_0529: 1.0, s_0794: 1.0}, 'YJR019C')
    r_o001 = add_reaction('r_o001', 'DHA synthesis', 0, 1000, 
                          {s_1260: -1.0, s_1101: -4.0, s_1212: -5.0, s_0794: -5.0, s_1275: -5.0, s_o001: 1.0, s_0529: 4.0, s_1203: 5.0, s_0803: 4.0}, 'DHAS')
    r_o002 = add_reaction('r_o002', 'DHA transport', -1000, 1000, {s_o001: -1.0, s_o002: 1.0})
    r_o003 = add_reaction('r_o003', 'DHA exchange', 0, 1000, {s_o002: -1.0})
    cobra_model.add_reactions([r_o001, r_o002, r_o003])
    return cobra_model

def create_gpr_model(model_file, env_file='evo_envs.csv', type='aromatic'):
    reframed_model_file = f'{os.path.dirname(model_file)}/reframed-GEM-{type}.xml'
    env_file = f'{os.path.dirname(model_file)}/evo_envs.csv'

    if not os.path.exists(reframed_model_file):
        cobra_model = cobra.io.read_sbml_model(model_file)
        # cobra_model.reactions.get_by_id('r_1714').lower_bound = -1000
        # cobra_model.reactions.get_by_id('r_1808').lower_bound = -1000
        if type == 'aromatic':
            cobra_model = add_aromatic_pathway(cobra_model)
        elif type == 'indigoidine':
            cobra_model = add_indigoidine_pathway(cobra_model)
        elif type == 'bikaverin':
            cobra_model = add_bikaverin_pathway(cobra_model)
        elif type == 'omega3':
            cobra_model = add_omega3_pathway(cobra_model)
        cobra.io.write_sbml_model(cobra_model, reframed_model_file)

    reframed_model = reframed.load_cbmodel(reframed_model_file, flavor='fbc2')
    # '''model_flux_basis_wt'''
    # reframed_model.reactions.R_r_1714.lb = -10  # D-glucose exchange
    # reframed_model.reactions.R_r_1654.lb = -1000  # ammonium exchange
    # reframed_model.reactions.R_r_2111.lb = 1  # growth
    # reframed_model.reactions.R_r_2111.ub = 1  # growth
    # '''mimic anaerobiosis'''
    # reframed_model.reactions.R_r_0438.ub = 0  # ferrocytochrome-c:oxygen oxidoreductase
    # reframed_model.reactions.R_r_1021.lb = 0  # succinate dehydrogenase (ubiquinone-6)
    # reframed_model.reactions.R_r_1021.ub = 0  # succinate dehydrogenase (ubiquinone-6)
    # reframed_model.reactions.R_r_0226.ub = 0  # ATP synthase
    # reframed_model.reactions.R_r_0662.ub = 0  # isocitrate lyase
    '''close threonine aldolase assumed to have minor role during growth in minimal medium,
        close to equilibrium -7.2 +- 4.4 [kJ/mol] (http://equilibrator.weizmann.ac.il)'''
    reframed_model.reactions.R_r_1040.ub = 0  # threonine aldolase
    '''close FMN reductases with relevance in apoptosis if these are on,
        anaerobic conditions generate in simulations succinate and no glycerol this does not match with the real biological state'''
    reframed_model.reactions.R_r_0441.ub = 0  # FMN reductase
    reframed_model.reactions.R_r_0442.ub = 0  # FMN reductase
    # '''close glycerol utilization pathway as glucose is used as c-source'''
    # reframed_model.reactions.R_r_0487.ub = 0  # glycerol dehydrogenase (NADP-dependent)
    # reframed_model.reactions.R_r_0488.ub = 0  # glycerol kinase
    '''close higher alcohol acetate ester -esterases'''
    reframed_model.reactions.R_r_0656.ub = 0  # isoamyl acetate-hydrolyzing esterase
    reframed_model.reactions.R_r_0657.ub = 0  # isobutyl acetate-hydrolyzing esterase
    '''close phenylacetaldehyde exchange assuming redox status preferring further conversion to phenyl ethanol synthesis'''
    reframed_model.reactions.R_r_2001.ub = 0  # phenylacetaldehyde exchange
    # '''constrain isocitrate dehydrogenase in NADPH generating direction'''
    # reframed_model.reactions.R_r_0659.lb = 0  # isocitrate dehydrogenase (NADP)
    # '''glycine synthesis from serine when glucose is the carbon source'''
    # reframed_model.reactions.R_r_0502.lb = 0  # glycine hydroxymethyltransferase
    # reframed_model.reactions.R_r_0503.lb = 0  # glycine hydroxymethyltransferase
    # reframed_model.reactions.R_r_0732.lb = 0  # methylenetetrahydrofolate dehydrogenase (NADP)
    # reframed_model.reactions.R_r_0733.lb = 0  # methylenetetrahydrofolate dehydrogenase (NADP)
    reframed_model.reactions.R_r_4581.ub = 0
    reframed_model.reactions.R_r_4582.ub = 0

    '''close all environmental uptake reactions but set them reversible (able to open for a specific environment later)'''
    env_df = pd.read_csv(env_file, index_col=0)
    uptake_reactions = set()
    for ind, row in env_df.iterrows():
        rxns = ast.literal_eval(env_df.loc[ind, 'reactions'])
        for rxn in rxns:
            uptake_reactions.add(rxn)
            reframed_model.reactions[rxn].lb = 0
            # reframed_model.reactions[rxn].ub = inf
            reframed_model.reactions[rxn].reversible = True

    gpr_model = reframed.core.transformation.gpr_transform(reframed_model, inplace=False, add_proteome=True, gene_prefix='G_', usage_prefix='u_', pseudo_genes=None)
    solver = reframed.solvers.solver_instance(gpr_model)
    solver.add_constraint('flux_ratio_constraint1', {'R_r_1589': 2, 'R_r_2000': -1})
    solver.add_constraint('flux_ratio_constraint2', {'R_r_1580': 2, 'R_r_1581': -1})
    solver.add_constraint('flux_ratio_constraint3', {'R_r_1865': 2, 'R_r_1862': -1})

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

    # u_desired_trait = set()
    # for r in desired_trait:
    #     for g in reframed_model.reactions[r].get_genes():
    #         # u_desired_trait.add('u_' + g[2:])
    #         print(r, reframed_model.reactions[r].name, 'u_' + g[2:], reframed_model.reactions[r].reversible)
    #     if not reframed_model.reactions[r].reversible:
    #         u_desired_trait.add(r)
    # u_desired_trait = list(u_desired_trait)
    # # for x in u_desired_trait:
    # #     print(x)

    return gpr_model, reframed_model


def filter_desired_trait(reframed_model, desired_trait):
    new_desired_trait = []
    for rxn in desired_trait:
        rxn_name = reframed_model.reactions[rxn].name
        if not (rxn_name.endswith('transport') or rxn_name.endswith('diffusion')):
            new_desired_trait.append(rxn)
    # print(new_desired_trait)
    return new_desired_trait
