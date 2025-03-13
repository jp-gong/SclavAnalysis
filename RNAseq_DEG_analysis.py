import os
import joblib
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

# Create mapping from SRR IDs to sample annotations
SRR_to_SampleAnnot = {}
rna_sample_file = os.path.join("RNAseqData", "RNAseqMetadata.xlsx")
RNAseqSample_df = pd.read_excel(rna_sample_file, index_col=0)
for sample_annotation, cur_df in RNAseqSample_df.groupby("SampleAnnotation"):
    cur_df = cur_df.sort_values(by="Run")
    assert cur_df.shape[0] == 2
    row1, row2 = cur_df.iloc[0], cur_df.iloc[1]
    annot = sample_annotation.replace(" ", "_")
    SRR_to_SampleAnnot[row1["Run"]] = annot + "_1"
    SRR_to_SampleAnnot[row2["Run"]] = annot + "_2"

# Load count data and preprocess
counts_file = os.path.join("RNAseqData", "total_expression.csv")
counts_df = pd.read_csv(counts_file, index_col=0)
counts_df = counts_df.rename(columns=SRR_to_SampleAnnot)
counts_df = counts_df.fillna(0).astype(int).T

# Create metadata DataFrame
metadata = pd.DataFrame(index=counts_df.index)
metadata["strain"] = metadata.index.str.split("_").str[0]
metadata["time"] = metadata.index.str.split("_").str[1]
metadata["sample"] = metadata["strain"] + "-" + metadata["time"]

# Set up DESeq2 analysis
inference = DefaultInference(n_cpus=16)
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design_factors="sample",  # design_factors=["time", "strain"] is equivalent to ~ time + strain in R
    refit_cooks=True,
    inference=inference
)
dds.deseq2()

# Save normalized counts
norm_counts_df = pd.DataFrame(
    data=dds.layers["normed_counts"],
    index=dds.obs.index,
    columns=dds.var.index
).T
joblib.dump(norm_counts_df, os.path.join("data", "NormalizedGeneExprValue.joblib"))

# Generate fold change dictionary for each experimental sample vs control
FoldChange_Dict = {}
for strain in ["C1", "OR", "NL"]:
    for time in ["08h", "14h", "36h", "72h"]:
        exp_sample = f"{strain}-{time}"
        control_sample = f"WT-{time}"
        stat_res = DeseqStats(dds, contrast=["sample", exp_sample, control_sample], inference=inference, quiet=True)
        stat_res.summary(quiet=True)
        FoldChange_Dict[exp_sample] = stat_res.results_df
joblib.dump(FoldChange_Dict, os.path.join("data", "FoldChange_Dict.joblib"))

# Function to get DEGs based on adjusted p-value and log2 fold-change thresholds
def GetDEGs(exp, control, padj_cutoff=0.005, log2fc_cutoff=1, UpDown=None):
    stat_res = DeseqStats(dds, contrast=["sample", exp, control], inference=inference, quiet=True)
    stat_res.summary(quiet=True)
    res = stat_res.results_df
    if UpDown == "Up":
        deg_condition = (res.baseMean > 10) & (res.padj < padj_cutoff) & (res.log2FoldChange > log2fc_cutoff)
    elif UpDown == "Down":
        deg_condition = (res.baseMean > 10) & (res.padj < padj_cutoff) & (res.log2FoldChange < -log2fc_cutoff)
    elif UpDown in [None, "None"]:
        deg_condition = (res.baseMean > 10) & (res.padj < padj_cutoff) & (abs(res.log2FoldChange) > log2fc_cutoff)
    else:
        raise ValueError("UpDown must be 'Up', 'Down', or None")
    return res[deg_condition]

# Get DEGs using a p-value cutoff of 0.005
pval = 0.005
DEGs = {"Up": {}, "Down": {}}
for updown in ["Up", "Down"]:
    for strain in ["C1", "OR", "NL"]:
        for time in ["08h", "14h", "36h", "72h"]:
            exp_sample = f"{strain}-{time}"
            control_sample = f"WT-{time}"
            res = GetDEGs(exp_sample, control_sample, padj_cutoff=pval, UpDown=updown)
            DEGs[updown][exp_sample] = res
joblib.dump(DEGs, os.path.join("data", "DEGs.joblib"))