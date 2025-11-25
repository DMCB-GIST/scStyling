"""
ZINB-based scStyleBulk Generator

This script transforms bulk RNA-seq data into single-cell-style expression profiles
by fitting Zero-Inflated Negative Binomial (ZINB) distributions to single-cell data
and mapping bulk expression values to the learned single-cell distribution.

Author: [Your Name]
Date: 2025
"""

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
from scipy.stats import zscore, norm
import matplotlib.pyplot as plt
import umap
from sklearn.preprocessing import StandardScaler

np.random.seed(42)


# =============================================================================
# Utility Functions
# =============================================================================

def wrap_with_tqdm(executor, fn, args_list, desc="Processing"):
    """
    Wrap parallel execution with a tqdm progress bar.
    
    Args:
        executor: ProcessPoolExecutor instance
        fn: Function to apply to each argument
        args_list: List of arguments to process
        desc: Description for progress bar
        
    Returns:
        Iterator with tqdm wrapper
    """
    return tqdm(executor.map(fn, args_list), total=len(args_list), desc=desc)


def neg_log_likelihood(params, data):
    """
    Calculate negative log-likelihood for ZINB distribution.
    
    Args:
        params: Tuple of (r, p, pi) parameters
        data: Observed expression data
        
    Returns:
        Negative log-likelihood value
    """
    r, p, pi = params
    likelihoods = (1 - pi) * nbinom.pmf(data, r, p) + (pi * (data == 0))
    likelihoods = np.maximum(likelihoods, 1e-10)
    return -np.sum(np.log(likelihoods))


def zinb_mean_variance(r, p, pi):
    """
    Compute mean and variance of ZINB distribution.
    
    Args:
        r: Size parameter
        p: Probability parameter
        pi: Zero-inflation parameter
        
    Returns:
        Tuple of (mean, variance)
    """
    mean = (1 - pi) * (r * (1 - p) / p)
    variance = mean * (1 + ((1 - p) / p) * r + pi * (r * (1 - p) / p))
    return mean, variance


def init_random_zinb():
    """
    Initialize ZINB parameters with random values.
    
    Returns:
        Tuple of (r, p, pi) initial parameters
    """
    return (
        np.random.uniform(0.5, 5),
        np.random.uniform(0.1, 0.3),
        np.random.uniform(0.7, 0.9)
    )


def process_gene(args):
    """
    Optimize ZINB parameters and generate samples for a single gene.
    
    Args:
        args: Tuple of (gene_name, expression_data, num_samples)
        
    Returns:
        Dictionary containing optimized parameters and generated samples
    """
    gene, data, num_samples = args

    r_init, p_init, pi_init = init_random_zinb()
    initial_guess = [r_init, p_init, pi_init]
    bounds = [(1e-5, None), (1e-5, 1 - 1e-5), (1e-5, 1 - 1e-5)]

    result = minimize(
        neg_log_likelihood, 
        initial_guess, 
        args=(data,), 
        bounds=bounds, 
        method="L-BFGS-B"
    )

    r_opt, p_opt, pi_opt = result.x
    mean, variance = zinb_mean_variance(r_opt, p_opt, pi_opt)

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
    

def generate_zinb_sample(r, p, pi, n_samples=1):
    """
    Generate samples from ZINB distribution using loop-based approach.
    
    Args:
        r: Size parameter
        p: Probability parameter
        pi: Zero-inflation parameter
        n_samples: Number of samples to generate
        
    Returns:
        Array of generated samples
    """
    data = np.zeros(n_samples)
    for i in range(n_samples):
        if np.random.random() < pi:
            data[i] = 0
        else:
            data[i] = nbinom.rvs(r, p)
    return data


def generate_zinb_sample_vec(r, p, pi, n_samples):
    """
    Generate samples from ZINB distribution using vectorized approach.
    
    Args:
        r: Size parameter
        p: Probability parameter
        pi: Zero-inflation parameter
        n_samples: Number of samples to generate
        
    Returns:
        Array of generated samples
    """
    mask = np.random.random(n_samples) >= pi
    data = np.zeros(n_samples)
    data[mask] = nbinom.rvs(r, p, size=np.sum(mask))
    return data



def create_sc_rate_df(ratios, counts):
    """
    Create cumulative distribution DataFrame for single-cell expression bins.
    
    Args:
        ratios: Probability ratios for each expression value
        counts: Value counts for expression levels
        
    Returns:
        DataFrame with start_rate, end_rate, and value_rate columns
    """
    cumulative = np.cumsum(ratios)
    start_rate = np.insert(cumulative[:-1], 0, 0)
    end_rate = cumulative
    return pd.DataFrame(
        {
            "start_rate": start_rate,
            "end_rate": end_rate,
            "value_rate": counts.index,
        }
    )
    


def plot_combined_umap(bulk_data, singlecell_data, scStylebulk_data, 
                       drug="Erlotinib.png", figsize=(12, 10)):
    """
    Generate a single UMAP plot combining bulk, single-cell, and scStyleBulk data.
    
    Args:
        bulk_data: DataFrame with bulk RNA-seq expression data (samples x genes)
        singlecell_data: DataFrame with single-cell RNA-seq data (cells x genes)
        scStylebulk_data: DataFrame with scStyleBulk data (samples x genes)
        drug: Path to save the output figure
        figsize: Figure size as tuple (width, height)
    
    Returns:
        umap_embedding: Combined UMAP coordinates for all data types
    """
    
    # Align genes across all datasets
    common_genes = list(set(bulk_data.columns) & 
                       set(singlecell_data.columns) & 
                       set(scStylebulk_data.columns))
    
    bulk_aligned = bulk_data[common_genes]
    sc_aligned = singlecell_data[common_genes]
    scStyle_aligned = scStylebulk_data[common_genes]
    
    print(f"Number of common genes: {len(common_genes)}")
    print(f"Bulk samples: {bulk_aligned.shape[0]}")
    print(f"Single-cell samples: {sc_aligned.shape[0]}")
    print(f"scStyleBulk samples: {scStyle_aligned.shape[0]}")
    
    # Create labels for each data type
    bulk_labels = ['Bulk'] * bulk_aligned.shape[0]
    sc_labels = ['Single-cell'] * sc_aligned.shape[0]
    scStyle_labels = ['scStyleBulk'] * scStyle_aligned.shape[0]
    
    # Combine all data
    combined_data = pd.concat([bulk_aligned, sc_aligned, scStyle_aligned], 
                              axis=0, ignore_index=True)
    combined_labels = bulk_labels + sc_labels + scStyle_labels
    
    # Standardize features
    scaler = StandardScaler()
    combined_scaled = scaler.fit_transform(combined_data)
    
    # Apply UMAP
    print("Computing UMAP embedding...")
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, 
                        random_state=42, metric='euclidean')
    umap_embedding = reducer.fit_transform(combined_scaled)
    
    # Create plot
    plt.figure(figsize=figsize)
    
    # Define colors and markers for each data type
    colors = {'Bulk': '#8da9cf', 'Single-cell': '#e7c583', 'scStyleBulk': '#c7dda6'}
    markers = {'Bulk': 's', 'Single-cell': 'o', 'scStyleBulk': '^'}
    sizes = {'Bulk': 80, 'Single-cell': 30, 'scStyleBulk': 80}
    
    # Plot each data type
    for data_type in ['Single-cell', 'Bulk', 'scStyleBulk']:
        mask = [label == data_type for label in combined_labels]
        plt.scatter(
            umap_embedding[mask, 0],
            umap_embedding[mask, 1],
            c=colors[data_type],
            marker=markers[data_type],
            s=sizes[data_type],
            alpha=0.5,
            label=data_type,
            edgecolors='black',
            linewidths=0.5
        )
    
    plt.xlabel('UMAP 1', fontsize=14)
    plt.ylabel('UMAP 2', fontsize=14)
    plt.title('Combined UMAP: Bulk vs Single-cell vs scStyleBulk', fontsize=16, fontweight='bold')
    plt.legend(loc='best', fontsize=12, frameon=True, shadow=True)
    plt.grid(False)
    plt.tight_layout()
    
    # Save figure
    plt.savefig(drug, dpi=300, bbox_inches='tight')
    print(f"UMAP plot saved to: {drug}")
    
    # Return embedding and labels for further analysis
    return umap_embedding, combined_labels


# =============================================================================
# Main Pipeline
# =============================================================================

def main(args):
    """
    Main execution pipeline:
    1. Load bulk and single-cell RNA-seq data
    2. Fit ZINB distributions to single-cell data
    3. Generate simulated single-cell data
    4. Transform bulk data into scStyleBulk format
    
    Args:
        args: Command-line arguments containing file paths and parameters
    """

    # -------------------------------------------------------------------------
    # Step 1: Load and preprocess data
    # -------------------------------------------------------------------------
    Bulk = pd.read_csv(f"{args.bulk_data}", index_col=0)
    SingleCell = pd.read_csv(f"{args.sc_data}", index_col=0)
    print("Load Bulk and Single-cell data")
    print('Bulk: ', Bulk.shape)
    print('Single-cell: ', SingleCell.shape)

    num_cores = max(1, os.cpu_count() - 2)
    print(f"Using {num_cores} cores for parallel processing")

    com_genes = list(set(Bulk.columns) & set(SingleCell.columns))
    Bulk_comgenes = Bulk[com_genes]
    SingleCell_comgenes = SingleCell[com_genes].round()
    num_samples = SingleCell_comgenes.shape[0]

    Bulk_comgenes.to_csv(f"{args.save_path}/Bulk_{args.drug}_comgenes.csv")
    SingleCell_comgenes.to_csv(f"{args.save_path}/SingleCell_{args.drug}_comgenes.csv")

    # -------------------------------------------------------------------------
    # Step 2: ZINB parameter estimation and sample generation
    # -------------------------------------------------------------------------
    print("Starting ZINB parameter estimation")
    
    gene_args = [(gene, SingleCell_comgenes[gene].values, num_samples) for gene in com_genes]

    all_results, all_samples = [], {}

    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        for result in wrap_with_tqdm(executor, process_gene, gene_args, desc="Processing"):
            if result is not None:
                all_results.append({k: v for k, v in result.items() if k != "samples"})
                all_samples[result["Gene"]] = result["samples"]

    simulated_sc_data = pd.DataFrame(all_samples)
    simulated_sc_data.to_csv(f"{args.save_path}/Simulated_singlecell_{args.drug}.csv")
    
    print("Completed ZINB parameter estimation")
    print(simulated_sc_data)
    
    print("Bulk shape:", Bulk.shape)
    print("Single-cell shape:", SingleCell.shape)
    print("Simulated Single-cell shape:", simulated_sc_data.shape)

    # -------------------------------------------------------------------------
    # Step 3: Transform bulk RNA-seq into scStyleBulk
    # -------------------------------------------------------------------------
    print("Starting transformation of bulk RNA-seq into scStyleBulk")
    
    z_Bulk = Bulk.apply(zscore, axis=0)
    cdf_Bulk = z_Bulk.apply(norm.cdf)

    all_genes = list(simulated_sc_data.keys())
    all_values = {gene: pd.Series(value) for gene, value in simulated_sc_data.items()}

    value_counts = {gene: values.value_counts().sort_index() for gene, values in all_values.items()}
    value_ratios = {gene: counts / counts.sum() for gene, counts in value_counts.items()}

    sc_rate_dfs = {gene: create_sc_rate_df(value_ratios[gene].values, value_counts[gene]) for gene in all_genes}

    bins = {gene: df["end_rate"].values for gene, df in sc_rate_dfs.items()}
    labels = {gene: df["value_rate"][1:].values for gene, df in sc_rate_dfs.items()}

    new_tmp_df = pd.DataFrame(
        {
            gene: pd.cut(
                cdf_Bulk[gene],
                bins=bins[gene],
                labels=labels[gene],
                include_lowest=True,
            )
            for gene in all_genes
        }
    )

    scStyleBulk = new_tmp_df.astype(float).fillna(0)

    output_path = f"{args.save_path}/scStyleBulk_{args.drug}.csv"
    scStyleBulk.to_csv(output_path)

    print("Completed Transformation")
    print("Saved scStyleBulk to:", output_path)
    print(scStyleBulk)
    
    
    plot_combined_umap(Bulk, SingleCell, scStyleBulk, drug=f"fig/{args.drug}.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Transform bulk RNA-seq data into single-cell-style expression profiles using ZINB modeling"
    )

    parser.add_argument("--drug", type=str, required=True, 
                        help="Drug name")
    parser.add_argument("--bulk_data", type=str, default="", 
                        help="Path to bulk RNA-seq data CSV file")
    parser.add_argument("--sc_data", type=str, default="", 
                        help="Path to single-cell RNA-seq data CSV file")
    parser.add_argument("--save_path", type=str, required=True, 
                        help="Directory path to save output files")

    args = parser.parse_args()
    main(args)