import os
import numpy as np
import pandas as pd
from pydeseq2.preprocessing import deseq2_norm
from utils import utils

# Create mapping from SRR to sample annotation
SRR_to_SampleAnnot = {}
rna_meta = pd.read_excel(os.path.join("RNAseqData", "RNAseqMetadata.xlsx"), index_col=0)
for sample_annot, group in rna_meta.groupby("SampleAnnotation"):
    group = group.sort_values(by="Run")
    assert group.shape[0] == 2, "Expected exactly two runs per sample"
    row1, row2 = group.iloc[0], group.iloc[1]
    annotation = sample_annot.replace(" ", "_")
    SRR_to_SampleAnnot[row1["Run"]] = f"{annotation}_1"
    SRR_to_SampleAnnot[row2["Run"]] = f"{annotation}_2"

# Load and preprocess counts data
counts_df = pd.read_csv(os.path.join("RNAseqData", "total_expression.csv"), index_col=0)
counts_df = counts_df.rename(columns=SRR_to_SampleAnnot).fillna(0).astype(int).T

# Normalize counts using DESeq2 normalization
normed_df, size_factor = deseq2_norm(counts_df)
normed_df = normed_df.T

# Reorder columns and apply log2 transformation
column_order = [
    'WT_08h', 'WT_14h', 'WT_36h', 'WT_72h',
    'C1_08h', 'C1_14h', 'C1_36h', 'C1_72h',
    'OR_08h', 'OR_14h', 'OR_36h', 'OR_72h',
    'NL_08h', 'NL_14h', 'NL_36h', 'NL_72h'
]
# Create column order with replicates (_1 and _2)
column_order = [f"{annot}_{i}" for annot in column_order for i in [1, 2]]
normed_df = normed_df[column_order]
normed_df = np.log2(normed_df + 1)

# Remove deleted genes from the normalized data
gene_df = pd.read_csv(os.path.join("Genomes", "Sclav_Gene_Info.csv"))
gene_df = gene_df.query("contig == 'plasmid'").query("(804919 > start) or (end > 1407036)")
deleted_genes = gene_df.locus_tag.tolist()
normed_df = normed_df.drop(deleted_genes, axis=0, errors='ignore')

# Perform ICA to obtain S and A matrices
S_final, A_final, stats_df = utils.perform_ICA(normed_df)

# Save the resulting matrices
S_final.to_csv("ICA_result_Matrix_S.csv")
A_final.to_csv("ICA_result_Matrix_A.csv")
