import os
# import sys
import time
from simulation import *
# from adage.robustness_analysis import *
from modeling import *
from concurrent.futures import ProcessPoolExecutor, wait
from config import *


WK_NO = min(10, os.cpu_count() - 2)


def run_dha01(output_dir):
    gpr_model, _ = create_gpr_model(model_xml, type='omega3')
    solution = simulate_engineered_strain(gpr_model, ynb_medium, [substrate_rxns['Glc']], 0.3, 10, ['R_r_o003'], evo_base_constr)
    solution.to_dataframe().to_csv(f'{output_dir}/dha01_results.csv')
    return solution.values


def run_dha02(ref_values, output_dir):
    gpr_model, _ = create_gpr_model(model_xml, type='omega3')
    solution = simulate_engineered_strain(gpr_model, ynb_medium, [substrate_rxns['Glc']], 0.3, 10, ['R_r_o003'], evo_base_constr | {'u_YGL055W': (2*ref_values['u_YGL055W'], 1000)})
    solution.to_dataframe().to_csv(f'{output_dir}/dha02_results.csv')


def run_dha03(ref_values, output_dir):
    gpr_model, _ = create_gpr_model(model_xml, type='omega3')
    solution = simulate_engineered_strain(gpr_model, ynb_medium, [substrate_rxns['Glc']], 0.3, 10, ['R_r_o003'], evo_base_constr | {'u_YOR245C': (0, 0)})
    solution.to_dataframe().to_csv(f'{output_dir}/dha03_results.csv')


def run_dha04(ref_values, output_dir):
    gpr_model, _ = create_gpr_model(model_xml, type='omega3')
    solution = simulate_engineered_strain(gpr_model, ynb_medium, [substrate_rxns['Glc']], 0.3, 10, ['R_r_o003'], evo_base_constr | {'u_YNR016C': (2*ref_values['u_YNR016C'], 1000)})
    solution.to_dataframe().to_csv(f'{output_dir}/dha04_results.csv')


def run_dha05(ref_values, output_dir):
    gpr_model, _ = create_gpr_model(model_xml, type='omega3')
    solution = simulate_engineered_strain(gpr_model, ynb_medium, [substrate_rxns['Glc']], 0.3, 10, ['R_r_o003'], evo_base_constr | {'u_YGL205W': (0, 0)})
    solution.to_dataframe().to_csv(f'{output_dir}/dha05_results.csv')


def run_dha06(ref_values, output_dir):
    gpr_model, _ = create_gpr_model(model_xml, type='omega3')
    solution = simulate_engineered_strain(gpr_model, ynb_medium, [substrate_rxns['Glc']], 0.3, 10, ['R_r_o003'], evo_base_constr | {'u_YKL182W': (2*ref_values['u_YKL182W'], 1000)})
    solution.to_dataframe().to_csv(f'{output_dir}/dha06_results.csv')


def run_dha07(ref_values, output_dir):
    gpr_model, _ = create_gpr_model(model_xml, type='omega3')
    solution = simulate_engineered_strain(gpr_model, ynb_medium, [substrate_rxns['Glc']], 0.3, 10, ['R_r_o003'], evo_base_constr | {'u_YNR008W': (0, 0)})
    solution.to_dataframe().to_csv(f'{output_dir}/dha07_results.csv')


if __name__ == '__main__':
    start_time = time.time()
    
    output_dir = 'output/omega3'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    ref_values = run_dha01(output_dir)
    futures = []
    with ProcessPoolExecutor(max_workers=WK_NO) as executor:
        futures.append(executor.submit(run_dha02, ref_values, output_dir))
        futures.append(executor.submit(run_dha03, ref_values, output_dir))
        futures.append(executor.submit(run_dha04, ref_values, output_dir))
        futures.append(executor.submit(run_dha05, ref_values, output_dir))
        futures.append(executor.submit(run_dha06, ref_values, output_dir))
        futures.append(executor.submit(run_dha07, ref_values, output_dir))
        wait(futures)
    # run_dha02(ref_values, output_dir)
    # run_dha03(ref_values, output_dir)
    # run_dha04(ref_values, output_dir)
    # run_dha05(ref_values, output_dir)
    # run_dha06(ref_values, output_dir)
    # run_dha07(ref_values, output_dir)

    print(f"Completed in {(time.time() - start_time) / 60 / 60:.2f}h.\n")