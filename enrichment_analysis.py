import matplotlib
matplotlib.use("Agg")  # Use non-interactive backend
import os
import pandas as pd
from collections import Counter
import joblib

from utils import GoAnalysis, CuratedKEGGAnalysis, utils

# ---------------------------------------------------------------------------
# Load DEGs joblib file
# ---------------------------------------------------------------------------
print("Loading DEGs data...")
degs_path = os.path.join("data", "DEGs.joblib")
DEGs = joblib.load(degs_path)
print("DEGs data loaded.\n")

# ---------------------------------------------------------------------------
# Step 1: GO Enrichment Analysis for Up-Regulated Genes
# ---------------------------------------------------------------------------
print("Step 1: Performing GO Enrichment Analysis for Up-Regulated Genes...")
go_results = {}
for strain in ["C1", "OR", "NL"]:
    for time in ["08h", "14h", "36h", "72h"]:
        exp_sample = f"{strain}-{time}"
        control_sample = f"WT-{time}"  # Defined but not used here
        res = DEGs["Up"][exp_sample]
        go_result = GoAnalysis.GO_Enrichment_Analysis(res)
        filter_condition = (
            (go_result["# of gene in GO"] >= 5)
            & (go_result["p value"] < 0.05)
            & (go_result["# of DEGs in GO"] >= 3)
            & (go_result["domain"] == "process")
            & (go_result["# of DEGs in GO"] / go_result["# of gene in GO"] > 200 / 6773)
        )
        go_results[exp_sample] = go_result[filter_condition].copy()

# Count GO term occurrences across samples
go_names = [df["GO name"].to_list() for df in go_results.values()]
go_name_counts = Counter(term for sublist in go_names for term in sublist)
duplicated_go_names = [name for name, count in go_name_counts.items() if count >= 2]

# Build DataFrames for p values and generatio (DEG ratio)
p_value = pd.DataFrame(index=duplicated_go_names)
generatio = pd.DataFrame(index=duplicated_go_names)
for sample, df in go_results.items():
    filtered_df = df[df["GO name"].isin(duplicated_go_names)][["p value", "GO name"]]
    filtered_df = filtered_df.set_index("GO name").rename(columns={"p value": sample})
    filtered_df = filtered_df[filtered_df[sample] < 0.05]
    p_value = pd.concat([p_value, filtered_df], axis=1)

    df["generatio"] = df["# of DEGs in GO"] / df["# of gene in GO"]
    filtered_df = df[df["GO name"].isin(duplicated_go_names)][["generatio", "GO name"]]
    filtered_df = filtered_df.set_index("GO name").rename(columns={"generatio": sample})
    generatio = pd.concat([generatio, filtered_df], axis=1)

# Clean index labels
p_value.index = p_value.index.str.title().str.replace("'", "")
generatio.index = generatio.index.str.title().str.replace("'", "")
p_value.sort_index(inplace=True)
generatio.sort_index(inplace=True)

selected_indexes = [
    "Translation",
    "Arginine Biosynthetic Process",
    "Amino Acid Transport",
    "Protein Folding",
    "Thiamine Biosynthetic Process",
]
p_value = p_value.loc[selected_indexes]
generatio = generatio.loc[selected_indexes]

# Draw the ontology plot for Up-regulated genes (figure saved only)
utils.DrawOntologyPlot(p_value, generatio, figsize=(5, 8),
                         save_path="images/GO_Enrichment_Analysis_Up.svg")
print("Step 1 completed.\n")

# ---------------------------------------------------------------------------
# Step 2: GO Enrichment Analysis for Down-Regulated Genes
# ---------------------------------------------------------------------------
print("Step 2: Performing GO Enrichment Analysis for Down-Regulated Genes...")
go_results = {}
for strain in ["C1", "OR", "NL"]:
    for time in ["08h", "14h", "36h", "72h"]:
        exp_sample = f"{strain}-{time}"
        control_sample = f"WT-{time}"
        res = DEGs["Down"][exp_sample]
        go_result = GoAnalysis.GO_Enrichment_Analysis(res)
        filter_condition = (
            (go_result["# of gene in GO"] >= 5)
            & (go_result["p value"] < 0.05)
            & (go_result["# of DEGs in GO"] >= 3)
            & (go_result["domain"] == "process")
            & (go_result["# of DEGs in GO"] / go_result["# of gene in GO"] > 150 / 6773)
        )
        go_results[exp_sample] = go_result[filter_condition].copy()

# Count GO term occurrences
go_names = [df["GO name"].to_list() for df in go_results.values()]
go_name_counts = Counter(term for sublist in go_names for term in sublist)
duplicated_go_names = [name for name, count in go_name_counts.items() if count >= 2]

p_value = pd.DataFrame(index=duplicated_go_names)
generatio = pd.DataFrame(index=duplicated_go_names)
for sample, df in go_results.items():
    filtered_df = df[df["GO name"].isin(duplicated_go_names)][["p value", "GO name"]]
    filtered_df = filtered_df.set_index("GO name").rename(columns={"p value": sample})
    filtered_df = filtered_df[filtered_df[sample] < 0.05]
    p_value = pd.concat([p_value, filtered_df], axis=1)

    df["generatio"] = df["# of DEGs in GO"] / df["# of gene in GO"]
    filtered_df = df[df["GO name"].isin(duplicated_go_names)][["generatio", "GO name"]]
    filtered_df = filtered_df.set_index("GO name").rename(columns={"generatio": sample})
    generatio = pd.concat([generatio, filtered_df], axis=1)

p_value.index = p_value.index.str.title().str.replace("'", "")
generatio.index = generatio.index.str.title().str.replace("'", "")
p_value.sort_index(inplace=True)
generatio.sort_index(inplace=True)

selected_indexes = [
    "Proton Motive Force-Driven Atp Synthesis",
    "Base-Excision Repair",
    "Translation",
    "Protein Folding",
    "Carbohydrate Metabolic Process",
]
p_value = p_value.loc[selected_indexes]
generatio = generatio.loc[selected_indexes]

utils.DrawOntologyPlot(p_value, generatio, figsize=(5, 8),
                         save_path="images/GO_Enrichment_Analysis_Down.svg")
print("Step 2 completed.\n")

# ---------------------------------------------------------------------------
# Step 3: KEGG Enrichment Analysis for Up-Regulated Genes
# ---------------------------------------------------------------------------
print("Step 3: Performing KEGG Enrichment Analysis for Up-Regulated Genes...")
kegg_up_results = {}
for strain in ["C1", "OR", "NL"]:
    for time in ["08h", "14h", "36h", "72h"]:
        exp_sample = f"{strain}-{time}"
        control_sample = f"WT-{time}"
        res = DEGs["Up"][exp_sample]
        kegg_result = CuratedKEGGAnalysis.KEGG_Enrichment_Analysis(res)
        filter_condition = (
            (kegg_result["# of gene in Pathway"] >= 5)
            & (kegg_result["p value"] < 0.05)
        )
        kegg_up_results[exp_sample] = kegg_result[filter_condition].copy()

pathway_names = [df["KEGGterm name"].to_list() for df in kegg_up_results.values()]
kegg_name_counts = Counter(term for sublist in pathway_names for term in sublist)
duplicated_kegg_names = [name for name, count in kegg_name_counts.items() if count >= 1]

p_value = pd.DataFrame(index=duplicated_kegg_names)
generatio = pd.DataFrame(index=duplicated_kegg_names)
for sample, df in kegg_up_results.items():
    filtered_df = df[df["KEGGterm name"].isin(duplicated_kegg_names)][["p value", "KEGGterm name"]]
    filtered_df = filtered_df.set_index("KEGGterm name").rename(columns={"p value": sample})
    filtered_df = filtered_df[filtered_df[sample] < 0.05]
    p_value = pd.concat([p_value, filtered_df], axis=1)

    df["generatio"] = df["# of DEGs in Pathways"] / df["# of gene in Pathway"]
    filtered_df = df[df["KEGGterm name"].isin(duplicated_kegg_names)][["generatio", "KEGGterm name"]]
    filtered_df = filtered_df.set_index("KEGGterm name").rename(columns={"generatio": sample})
    generatio = pd.concat([generatio, filtered_df], axis=1)

p_value.index = p_value.index.str.title().str.replace("'", "")
generatio.index = generatio.index.str.title().str.replace("'", "")
p_value.sort_index(inplace=True)
generatio.sort_index(inplace=True)

selected_indexes = [
    "Clavulanic Acid Biosynthesis",
    "Cephamycin C & Cephalosporin C Biosynthesis",
    "Arginine Biosynthesis",
    "Ornithine Biosynthesis, Glutamate => Ornithine",
    "Ribosome",
    "Assimilatory Sulfate Reduction, Sulfate => H2S",
    "Sulfur Metabolism",
    "Quorum Sensing",
    "Pyruvate Metabolism",
    "Leucine Biosynthesis, 2-Oxoisovalerate => 2-Oxoisocaproate",
    "Valine, Leucine And Isoleucine Biosynthesis",
]
p_value = p_value.loc[selected_indexes]
generatio = generatio.loc[selected_indexes]

utils.DrawOntologyPlot(p_value, generatio, figsize=(5, 8),
                         save_path="images/KEGG_Enrichment_Analysis_Up.svg")
print("Step 3 completed.\n")

# ---------------------------------------------------------------------------
# Step 4: KEGG Enrichment Analysis for Down-Regulated Genes
# ---------------------------------------------------------------------------
print("Step 4: Performing KEGG Enrichment Analysis for Down-Regulated Genes...")
kegg_down_results = {}
for strain in ["C1", "OR", "NL"]:
    for time in ["08h", "14h", "36h", "72h"]:
        exp_sample = f"{strain}-{time}"
        control_sample = f"WT-{time}"
        res = DEGs["Down"][exp_sample]
        kegg_result = CuratedKEGGAnalysis.KEGG_Enrichment_Analysis(res)
        filter_condition = (
            (kegg_result["# of gene in Pathway"] >= 5)
            & (kegg_result["p value"] < 0.05)
        )
        kegg_down_results[exp_sample] = kegg_result[filter_condition].copy()

pathway_names = [df["KEGGterm name"].to_list() for df in kegg_down_results.values()]
kegg_name_counts = Counter(term for sublist in pathway_names for term in sublist)
duplicated_kegg_names = [name for name, count in kegg_name_counts.items() if count >= 1]

p_value = pd.DataFrame(index=duplicated_kegg_names)
generatio = pd.DataFrame(index=duplicated_kegg_names)
for sample, df in kegg_down_results.items():
    filtered_df = df[df["KEGGterm name"].isin(duplicated_kegg_names)][["p value", "KEGGterm name"]]
    filtered_df = filtered_df.set_index("KEGGterm name").rename(columns={"p value": sample})
    filtered_df = filtered_df[filtered_df[sample] < 0.05]
    p_value = pd.concat([p_value, filtered_df], axis=1)

    df["generatio"] = df["# of DEGs in Pathways"] / df["# of gene in Pathway"]
    filtered_df = df[df["KEGGterm name"].isin(duplicated_kegg_names)][["generatio", "KEGGterm name"]]
    filtered_df = filtered_df.set_index("KEGGterm name").rename(columns={"generatio": sample})
    generatio = pd.concat([generatio, filtered_df], axis=1)

p_value.index = p_value.index.str.title().str.replace("'", "")
generatio.index = generatio.index.str.title().str.replace("'", "")
p_value.sort_index(inplace=True)
generatio.sort_index(inplace=True)

selected_indexes = [
    "Ribosome",
    "Pyrimidine Metabolism",
    "Citrate Cycle (Tca Cycle)",
    "Cysteine And Methionine Metabolism",
    "Glycolysis / Gluconeogenesis",
    "Pentose Phosphate Pathway",
    "Purine Metabolism",
    "Ubiquinone And Other Terpenoid-Quinone Biosynthesis",
    "Staurosporine Biosynthesis",
    "Holomycin Bgc",
    "Naringenin Bgc",
]
p_value = p_value.loc[selected_indexes]
generatio = generatio.loc[selected_indexes]

utils.DrawOntologyPlot(p_value, generatio, figsize=(5, 8),
                         save_path="images/KEGG_Enrichment_Analysis_Down.svg")
print("Step 4 completed.\n")

print("All steps completed successfully!")
# ---------------------------------------------------------------------------