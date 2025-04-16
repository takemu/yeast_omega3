from modeling import *
from reframed.solvers.solver import VarType


# env_base = []
# proxy_fitness = ['R_r_1714', 'R_r_1654']  # reframed_rxns
# growth_weight = 1
# min_growth = 1e-6  # h^-1
# target_growth = 0.4  # h^-1
# combined = True


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
        for gpr_rxn in gpr_reactions([rxn], excludes=['_f', '_b']):
            gpr_constrs.update({gpr_rxn: bounds})
        for gpr_rxn in gpr_reactions([rxn], includes=['_f']):
            gpr_constrs.update({gpr_rxn: (0, bounds[1])})
        for gpr_rxn in gpr_reactions([rxn], includes=['_b']):
            gpr_constrs.update({gpr_rxn: (0, -bounds[0])})
    return gpr_constrs


def moma_solver_instance(gpr_model, reference, alg='moma'):
    solver = reframed.solvers.solver_instance(gpr_model)
    if alg == 'lmoma':
        if not hasattr(solver, 'lMOMA_flag'):
            solver.lMOMA_flag = True
            for r_id in reference.keys():
                d_pos, d_neg = r_id + '_d+', r_id + '_d-'
                solver.add_variable(d_pos, 0, inf, update=False)
                solver.add_variable(d_neg, 0, inf, update=False)
            solver.update()
            for r_id in reference.keys():
                d_pos, d_neg = r_id + '_d+', r_id + '_d-'
                solver.add_constraint('c' + d_pos, {r_id: -1, d_pos: 1}, '>', -reference[r_id], update=False)
                solver.add_constraint('c' + d_neg, {r_id: 1, d_neg: 1}, '>', reference[r_id], update=False)
            solver.update()
    elif alg == 'room':
        U = 1e6
        L = -1e6
        delta = 0.1
        epsilon = 0.001
        if not hasattr(solver, 'ROOM_flag'):
            solver.ROOM_flag = True

            for r_id in reference.keys():
                y_i = 'y_' + r_id
                solver.add_variable(y_i, 0, 1, vartype=VarType.BINARY, update=False)
            solver.update()

            for r_id in reference.keys():
                y_i = 'y_' + r_id
                if isinstance(reference[r_id], tuple) or isinstance(reference[r_id], list):
                    w_i_min = reference[r_id][0] if reference[r_id][0] != -inf else -1000
                    w_i_max = reference[r_id][1] if reference[r_id][1] != inf else 1000
                else:
                    w_i_min = reference[r_id]
                    w_i_max = reference[r_id]
                w_u = w_i_max + delta * abs(w_i_max) + epsilon
                w_l = w_i_min - delta * abs(w_i_min) - epsilon
                solver.add_constraint('c' + r_id + '_u', {r_id: 1, y_i: (w_u - U)}, '<', w_u, update=False)
                solver.add_constraint('c' + r_id + '_l', {r_id: 1, y_i: (w_l - L)}, '>', w_l, update=False)
            solver.update()

    return solver


def compare_values(value1, value2):
    if value2 - value1 > 0:
        return 'up'
    elif value2 - value1 < 0:
        return 'down'
    else:
        return 'unchanged'


def sum_flux(fluxes, gpr_rxns, abs_value=False, gpr_signs=None):
    # if abs_value:
    #     # return fluxes[gpr_rxns[0]].abs().sum() + np.sum(np.abs(fluxes[gpr_rxns[1]].values - fluxes[gpr_rxns[2]].values))
    #     return fluxes[gpr_rxns[0]].abs().sum() + np.sum(np.maximum(gpr_signs * (fluxes[gpr_rxns[1]].values - fluxes[gpr_rxns[2]].values), 0))
    # else:
    return fluxes[gpr_rxns[0]].sum() + fluxes[gpr_rxns[1]].sum() - fluxes[gpr_rxns[2]].sum()