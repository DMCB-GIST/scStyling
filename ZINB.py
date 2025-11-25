import os
import random
import argparse
from functools import partial
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm
from scipy.stats import nbinom
from scipy.optimize import minimize

# Fix random seed for reproducibility
np.random.seed(42)


# ------------------------------
# Utility Functions
# ------------------------------

def wrap_with_tqdm(executor, fn, args_list, desc="Processing"):
    """Wrap parallel execution with tqdm progress bar"""
    return tqdm(executor.map(fn, args_list), total=len(args_list), desc=desc)


def neg_log_likelihood(params, data):
    """Calculate negative log-likelihood for ZINB distribution"""
    r, p, pi = params
    likelihoods = (1 - pi) * nbinom.pmf(data, r, p) + (pi * (data == 0))
    likelihoods = np.maximum(likelihoods, 1e-10)  # avoid log(0)
    return -np.sum(np.log(likelihoods))


def zinb_mean_variance(r, p, pi):
    """Compute mean and variance of ZINB distribution"""
    mean = (1 - pi) * (r * (1 - p) / p)
    variance = mean * (1 + ((1 - p) / p) * r + pi * (r * (1 - p) / p))
    return mean, variance


def init_random_zinb():
    """Initialize ZINB parameters with random values"""
    return (
        np.random.uniform(0.5, 5),    # r
        np.random.uniform(0.1, 0.3),  # p
        np.random.uniform(0.7, 0.9)   # pi
    )


def generate_zinb_sample(r, p, pi, n_samples=1):
    """Generate samples from ZINB distribution (loop-based)"""
    data = np.zeros(n_samples)
    for i in range(n_samples):
        if np.random.random() < pi:
            data[i] = 0
        else:
            data[i] = nbinom.rvs(r, p)
    return data


def generate_zinb_sample_vec(r, p, pi, n_samples):
    """Generate samples from ZINB distribution (vectorized version)"""
    mask = np.random.random(n_samples) >= pi
    data = np.zeros(n_samples)
    data[mask] = nbinom.rvs(r, p, size=np.sum(mask))
    return data


def process_gene(args):
    """Optimize ZINB parameters and generate samples for a single gene"""
    gene, data, num_samples = args

    # Initialize parameters
    r_init, p_init, pi_init = init_random_zinb()
    initial_guess = [r_init, p_init, pi_init]
    bounds = [(1e-5, None), (1e-5, 1 - 1e-5), (1e-5, 1 - 1e-5)]

    # Parameter optimization
    result = minimize(
        neg_log_likelihood, 
        initial_guess, 
        args=(data,), 
        bounds=bounds, 
        method="L-BFGS-B"
    )

    r_opt, p_opt, pi_opt = result.x
    mean, variance = zinb_mean_variance(r_opt, p_opt, pi_opt)

    # Generate samples with optimized parameters
    samples = generate_zinb_sample(r_opt, p_opt, pi_opt, num_samples)

    return {
        "Gene": gene,
        "Optimal r": r_opt,
        "Optimal p": p_opt,
        "Optimal pi": pi_opt,
        "Optimal mean": mean,
        "Optimal variance": variance,
        "samples": samples
    }


# ------------------------------
# Main Function
# ------------------------------

def main(args):
    """Main execution: load bulk and single-cell data, fit ZINB, generate samples"""

    # Load bulk and single-cell data
    Bulk = pd.read_csv(f"{args.bulk_data}",index_col=0)
    SingleCell = pd.read_csv(f"{args.sc_data}",index_col=0)
    print("Load Bulk and Single-cell data")
    print('Bulk: ', Bulk.shape)
    print('Single-cell: ', SingleCell.shape)

    # Determine available CPU cores for parallel processing
    num_cores = max(1, os.cpu_count() - 2)
    print(f"Using {num_cores} cores for parallel processing")

    # Extract common genes between bulk and single-cell data
    com_genes = list(set(Bulk.columns) & set(SingleCell.columns))
    Bulk_comgenes = Bulk[com_genes]
    SingleCell_comgenes = SingleCell[com_genes].round()
    num_samples = SingleCell_comgenes.shape[0]

    # Save preprocessed data
    Bulk_comgenes.to_csv(f"{args.save_path}/Bulk_{args.drug}.csv")
    SingleCell_comgenes.to_csv(f"{args.save_path}/SingleCell_{args.drug}.csv")

    # Prepare arguments for parallel processing
    gene_args = [(gene, SingleCell_comgenes[gene].values, num_samples) for gene in com_genes]

    all_results, all_samples = [], {}

    # Parallel processing
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        for result in wrap_with_tqdm(executor, process_gene, gene_args, desc="Processing"):
            if result is not None:
                all_results.append({k: v for k, v in result.items() if k != "samples"})
                all_samples[result["Gene"]] = result["samples"]

    # Save simulated single-cell data
    samples_df = pd.DataFrame(all_samples)
    samples_df.to_csv(f"{args.save_path}/Simulated_SingleCell_{args.drug}.csv")

    print(samples_df)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Estimate ZINB parameters")

    parser.add_argument("--drug", type=str, required=True, help="Prefix for output file names")
    parser.add_argument("--bulk_data", type=str, default="", help="Path to bulk RNA-seq data")
    parser.add_argument("--sc_data", type=str, default="", help="Path to single-cell RNA-seq data")
    parser.add_argument("--save_path", type=str, required=True, help="Path to save results")

    args = parser.parse_args()
    main(args)
