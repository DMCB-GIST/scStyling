import os
import argparse
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore, norm


def mean(SingleCell, generated_data, filename):
    """
    Compare gene-wise mean expression between original and generated single-cell data.
    Save a scatter plot showing the comparison.
    """
    # Calculate mean expression for each gene
    original_means = SingleCell.mean()
    generated_means = generated_data.mean()

    print("Original gene-wise means:\n", original_means)
    print("Generated gene-wise means:\n", generated_means)

    # Create scatter plot
    plt.figure(figsize=(10, 8))
    plt.scatter(
        original_means,
        generated_means,
        alpha=0.5,
        color="blue",
        s=20,
    )

    # Add y = x reference line
    max_val = max(original_means.max(), generated_means.max())
    plt.plot([0, max_val], [0, max_val], "r--", alpha=0.8)

    plt.title("Gene-wise Mean Expression Comparison", fontsize=12)
    plt.xlabel("Original Mean Expression", fontsize=10)
    plt.ylabel("Generated Mean Expression", fontsize=10)

    plt.grid(False)
    plt.tight_layout()

    # Save figure
    folder = "fig"
    os.makedirs(folder, exist_ok=True)
    new_filename = f"{folder}/{filename}.png"
    plt.savefig(new_filename, dpi=300, bbox_inches="tight")



def create_sc_rate_df(ratios, counts):
    """Create cumulative distribution DataFrame for a gene"""
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
    
        
def main(args):
    """
    Load bulk, single-cell, and simulated single-cell data.
    Perform distribution alignment (scStyleBulk transformation)
    and save the transformed bulk dataset.
    """

    # ---- Load data ----
    SingleCell = pd.read_csv(args.sc_data, index_col=0)
    Bulk = pd.read_csv(args.bulk_data, index_col=0)
    simulated_singlecell = pd.read_csv(args.simulated_sc_data, index_col=0)

    # Use common genes across all datasets
    com_genes = list(set(SingleCell.columns) & set(simulated_singlecell.columns) & set(Bulk.columns))
    SingleCell = SingleCell[com_genes]
    simulated_singlecell = simulated_singlecell[com_genes]
    Bulk = Bulk[com_genes]

    print("Loaded Bulk and Single-cell data")
    print("Bulk shape:", Bulk.shape)
    print("Single-cell shape:", SingleCell.shape)
    print("Simulated Single-cell shape:", simulated_singlecell.shape)

    # ---- Compare gene-wise mean ----
    mean(SingleCell, simulated_singlecell, filename=args.drug)

    # ---- Z-score normalization for Bulk ----
    z_Bulk = Bulk.apply(zscore, axis=0)
    cdf_Bulk = z_Bulk.apply(norm.cdf)

    # ---- Construct single-cell distribution ----
    all_genes = list(simulated_singlecell.keys())
    all_values = {gene: pd.Series(value) for gene, value in simulated_singlecell.items()}

    # Calculate probability distribution for each gene
    value_counts = {gene: values.value_counts().sort_index() for gene, values in all_values.items()}
    value_ratios = {gene: counts / counts.sum() for gene, counts in value_counts.items()}

    sc_rate_dfs = {gene: create_sc_rate_df(value_ratios[gene].values, value_counts[gene])for gene in all_genes}

    # ---- Map bulk CDF values to single-cell bins ----
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

    # Convert transformed bulk into numeric matrix
    scStyleBulk_Bulk = new_tmp_df.astype(float).fillna(0)

    # Save results
    output_path = f"{args.save_path}/scStyleBulk_{args.drug}.csv"
    scStyleBulk_Bulk.to_csv(output_path)
    print("Saved scStyleBulk-transformed bulk to:", output_path)
    print(scStyleBulk_Bulk)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform scStyleBulk transformation")
    parser.add_argument("--drug", type=str, required=True, help="Prefix for output file names")
    parser.add_argument("--bulk_data", type=str, default="", help="Path to bulk RNA-seq data")
    parser.add_argument("--sc_data", type=str, default="", help="Path to single-cell RNA-seq data")
    parser.add_argument("--simulated_sc_data", type=str, default="", help="Path to simulated single-cell RNA-seq data")
    parser.add_argument("--save_path", type=str, required=True, help="Directory to save results")

    args = parser.parse_args()
    main(args)
