from reframed.solvers.solver import VarType
from modeling import *
import adage.utils as utils
from adage.utils import *
import time
import numpy as np
import scipy as sp
from scipy import stats
from scipy.stats import *
import matplotlib.pyplot as plt
from reframed.solvers.solution import Status
from config import *
import itertools
import random


def mutation_sampling(single_mutations, n_mutations=5, n_samples=100):
    mutation_samples = set()
    i = 0
    while len(mutation_samples) < n_samples and i < 3000:
        comb = tuple(sorted(random.sample(single_mutations, min(n_mutations, len(single_mutations)))))
        mutation_samples.add(comb)
        i += 1
    return list(mutation_samples)


def robustness_analysis(gpr_model, fva_df, nutrients, desired_trait, n_mutations=5, n_samples=100):
    single_mutations = fva_df[(fva_df.index.str.startswith('u_')) & (fva_df['a_flux']!=fva_df['p_flux'])].index.tolist()
    # single_mutations = fva_df[(fva_df.index.str.startswith('y_')) & (fva_df['a_flux']>=0.5)].index
    # single_mutations = [rxn.replace('y_u_', 'u_') for rxn in single_mutations]

    dt_gpr_rxns = (gpr_reactions(desired_trait, excludes=['_f', '_b']), gpr_reactions(desired_trait, includes=['_f']), gpr_reactions(desired_trait, includes=['_b']))
    ut_gpr_rxns = (gpr_reactions(nutrients, excludes=['_f', '_b']), gpr_reactions(nutrients, includes=['_f']), gpr_reactions(nutrients, includes=['_b']))
    fva_df.loc['proxy_fitness', ['p_flux', 'a_flux']] = [fva_df.loc[r_biomass, 'p_flux'] / -sum_flux(fva_df['p_flux'], ut_gpr_rxns),
                                                         fva_df.loc[r_biomass, 'a_flux'] / -sum_flux(fva_df['a_flux'], ut_gpr_rxns)]
    fva_df.loc['desired_trait', ['p_flux', 'a_flux']] = [sum_flux(fva_df['p_flux'], dt_gpr_rxns) / -sum_flux(fva_df['p_flux'], ut_gpr_rxns), 
                                                         sum_flux(fva_df['a_flux'], dt_gpr_rxns) / -sum_flux(fva_df['a_flux'], ut_gpr_rxns)]

    ref_fluxes = fva_df.loc[gpr_model.u_reactions, 'p_flux'].to_dict()
    quad_obj = {(r_id, r_id): 1 for r_id in ref_fluxes.keys()}
    lin_obj = {r_id: -2 * ref_fluxes[r_id] for r_id in ref_fluxes.keys()}
    solver = reframed.solvers.solver_instance(gpr_model)
    results = []
    pro_fitness_single_mutations = []
    constraints = gpr_conversion({rxn: (-MAX_FLUX, 0) for rxn in nutrients}) | {r_biomass: (min_growth, MAX_FLUX)}
    for mutation in single_mutations:
        solution = solver.solve(lin_obj, 
                            quadratic=quad_obj, 
                            minimize=True, 
                            constraints=constraints | {mutation: (fva_df.loc[mutation, 'a_flux'], fva_df.loc[mutation, 'a_flux'])})
        if solution.status == Status.OPTIMAL or solution.status == Status.SUBOPTIMAL:
            pf = solution.values[r_biomass] / -sum_flux(solution.to_dataframe()['value'], ut_gpr_rxns)
            dt = sum_flux(solution.to_dataframe()['value'], dt_gpr_rxns) / -sum_flux(solution.to_dataframe()['value'], ut_gpr_rxns)
            if (pf - fva_df.loc['proxy_fitness', 'p_flux'] > 0) and (fva_df.loc['desired_trait', 'p_flux'] - dt >= 0):
                results.append([mutation, 1, pf, dt])
            # if (pf - fva_df.loc['proxy_fitness', 'p_flux']) / (fva_df.loc['proxy_fitness', 'a_flux'] - fva_df.loc['proxy_fitness', 'p_flux']) >= 0.05:
                pro_fitness_single_mutations.append(mutation)
    multiple_mutations = []
    for i in range(1, n_mutations):
        multiple_mutations += mutation_sampling(pro_fitness_single_mutations, i+1, n_samples)
    for mutation_comb in multiple_mutations:
        solution = solver.solve(lin_obj, 
                            quadratic=quad_obj, 
                            minimize=True, 
                            constraints=constraints | {rxn: (fva_df.loc[rxn, 'a_flux'], fva_df.loc[rxn, 'a_flux']) for rxn in mutation_comb})
        if solution.status == Status.OPTIMAL or solution.status == Status.SUBOPTIMAL:
            pf = solution.values[r_biomass] / -sum_flux(solution.to_dataframe()['value'], ut_gpr_rxns)
            dt = sum_flux(solution.to_dataframe()['value'], dt_gpr_rxns) / -sum_flux(solution.to_dataframe()['value'], ut_gpr_rxns)
            if (pf - fva_df.loc['proxy_fitness', 'p_flux'] > 0) and (fva_df.loc['desired_trait', 'p_flux'] - dt >= 0):
                results.append([mutation_comb, len(mutation_comb), pf, dt])
    
    ra_df = pd.DataFrame(results, columns=['mutations', 'n', 'proxy_fitness', 'desired_trait']).set_index('mutations')
    ra_df.loc['p_all', ['n', 'proxy_fitness', 'desired_trait']] = [0, fva_df.loc['proxy_fitness', 'p_flux'], fva_df.loc['desired_trait', 'p_flux']]
    ra_df.loc['a_all', ['n', 'proxy_fitness', 'desired_trait']] = [len(single_mutations), fva_df.loc['proxy_fitness', 'a_flux'], fva_df.loc['desired_trait', 'a_flux']]
    # ra_df = ra_df.rename(index=lambda x: str(x).replace('__45__', '-'))
    return ra_df