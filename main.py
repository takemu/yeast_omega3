import os
import sys
from adaptx.simulation import *
from adaptx.robustness_analysis import *
from modeling import *
from concurrent.futures import ProcessPoolExecutor, wait
from config import *


WK_NO = min(10, os.cpu_count() - 2)


def run_indigoidine(n_mutations, n_samples, output_dir):
    gpr_model, _ = create_gpr_model(model_xml, type='indigoidine')
    if not os.path.exists(f'{output_dir}/indigoidine_fva_results.csv'):
        r_fluxes = simulate_engineered_strain(gpr_model, ynb_medium, [substrate_rxns['Gal']], 0.2, 4, ['R_r_i003'], evo_base_constr)
        fva_df = simulate_adaptation(gpr_model, r_fluxes, [substrate_rxns['Gal']] + yp_medium)
        fva_df.to_csv(f'{output_dir}/indigoidine_fva_results.csv')
    else:
        fva_df = pd.read_csv(f'{output_dir}/indigoidine_fva_results.csv', index_col=0)
    ra_df = robustness_analysis(gpr_model, fva_df, [substrate_rxns['Gal']] + yp_medium, ['R_r_i003'], n_mutations, n_samples)
    ra_df.to_csv(f'{output_dir}/indigoidine_ra_results.csv')


def run_bikaverin(n_mutations, n_samples, output_dir):
    gpr_model, _ = create_gpr_model(model_xml, type='bikaverin')
    if not os.path.exists(f'{output_dir}/bikaverin_fva_results.csv'):
        r_fluxes = simulate_engineered_strain(gpr_model, ynb_medium, [substrate_rxns['Gal']], 0.15, 4, ['R_r_b003'], evo_base_constr)
        fva_df = simulate_adaptation(gpr_model, r_fluxes, [substrate_rxns['Gal']] + yp_medium)
        fva_df.to_csv(f'{output_dir}/bikaverin_fva_results.csv')
    else:
        fva_df = pd.read_csv(f'{output_dir}/bikaverin_fva_results.csv', index_col=0)
    ra_df = robustness_analysis(gpr_model, fva_df, [substrate_rxns['Gal']] + yp_medium, ['R_r_b003'], n_mutations, n_samples)
    ra_df.to_csv(f'{output_dir}/bikaverin_ra_results.csv')


if __name__ == '__main__':
    start_time = time.time()

    n_samples = 0
    for i, arg in enumerate(sys.argv):
        if i == 1:
            n_mutations = int(arg)
        elif i == 2:
            n_samples = int(arg)
    
    output_dir = 'output/omega3'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    futures = []
    with ProcessPoolExecutor(max_workers=WK_NO) as executor:
        futures.append(executor.submit(run_indigoidine, n_mutations, n_samples, output_dir))
        futures.append(executor.submit(run_bikaverin, n_mutations, n_samples, output_dir))
        wait(futures)

    print(f"Completed in {(time.time() - start_time) / 60 / 60:.2f}h.\n")