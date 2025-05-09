from modeling import *
from config import *
import numpy as np
from reframed.solvers import solver_instance
from reframed.cobra.simulation import FBA, pFBA
from math import floor, ceil


def gpr_reactions(reframed_rxns, includes=[], excludes=[]):
    u_reframed_rxns = []
    r_reframed_rxns = []
    for rxn in reframed_rxns:
        if rxn.startswith('u_'):
            u_reframed_rxns.append(rxn)
        elif rxn.startswith('R_'):
            r_reframed_rxns.append(rxn)
    gpr_r_rxns = sum([reframed_to_gpr_rxns[rxn] for rxn in r_reframed_rxns], [])
    for w in includes:
        gpr_r_rxns = [rxn for rxn in gpr_r_rxns if re.match(rf"R_r_.*{w}.*$", rxn)]
    for w in excludes:
        gpr_r_rxns = [rxn for rxn in gpr_r_rxns if not re.match(rf"R_r_.*{w}.*$", rxn)]
    return sorted(gpr_r_rxns + u_reframed_rxns)


def gpr_conversion(constraints):
    gpr_constrs = {}
    for rxn, bounds in constraints.items():
        gpr_rxns = gpr_reactions([rxn], excludes=['_f', '_b'])
        if gpr_rxns:
            for gpr_rxn in gpr_rxns:
                gpr_constrs.update({gpr_rxn: bounds})
        else:
            for gpr_rxn in gpr_reactions([rxn], includes=['_f']):
                gpr_constrs.update({gpr_rxn: (0, bounds[1])})
            for gpr_rxn in gpr_reactions([rxn], includes=['_b']):
                gpr_constrs.update({gpr_rxn: (0, -bounds[0])})
    return gpr_constrs


def simulate_engineered_strain(gpr_model, medium, carbon_source, target_growth, max_uptake, products, extra_constraints={}, ref_fluxes=None):
    """
    
    Arguments:
        gpr_model (CBModel): GPR transformed reframed cobra model
        medium (list): list of nutrients in medium
        carbon_source (list): list of substrates as carbon source
        target_growth (double): targeted specific growth rate (biomass)
        max_uptake (double): max uptake of carbon source
        products (list): list of target products of engineered strain
        extra_constraints (dict): extra constraints
        ref_fluxes (dict): fluxes of parental strain
    Returns:
        solution (Solution): fluxes of engineered strain

    """
    constraints = {rxn: (-MAX_FLUX, 0) for rxn in medium}
    constraints.update({rxn: (-max_uptake, 0) for rxn in carbon_source})
    constraints.update({r_biomass: (target_growth, MAX_FLUX)})
    constraints.update(extra_constraints)
    gpr_constraints = gpr_conversion(constraints)
    obj_rxns = gpr_reactions(products, includes=['_f']) + gpr_reactions(products, excludes=['_f', '_b'])
    solution = reframed.cobra.simulation.pFBA(gpr_model, 
                                                objective={rxn: 1 for rxn in obj_rxns}, 
                                                minimize=False, 
                                                constraints=gpr_constraints)
    return solution
