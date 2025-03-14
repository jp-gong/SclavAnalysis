import pandas as pd
import scipy.stats
import os
import re

def _ConstructKEGGPathways():
    with open(  "KEGG/ReconstructedPathway.txt" , "r" )  as f:
        lines = [line.strip() for line in f.read().split("\n")]

    to_be_update_pathway = {}
    pathways = list()
    KOGatheringState = False
    for line in lines:
        
        # Add Gene information
        if line.startswith("CRV15_RS"):
            to_be_update_pathway["gene_list"].extend( line.split(", ") )
            continue
        
        # Skip KO information
        if re.match('K[0-9]{5,5}', line):
            continue
        
        if KOGatheringState:
            to_be_update_pathway["# of gene in Pathway"] = len(to_be_update_pathway["gene_list"])
            pathways.append( to_be_update_pathway.copy() )
            KOGatheringState = False

        if line.startswith("##"):   to_be_update_pathway["Subunit"] = line[2:].strip()
        elif line.startswith("#"):     to_be_update_pathway["Unit"] = line[1:].strip()
        
        # Pathway information
        elif re.match('[0-9]{5,5}', line):
            to_be_update_pathway["Pathway code"] = line[:5]
            to_be_update_pathway["Pathway name"] = line[6:].split("(")[0].strip()
            to_be_update_pathway["# of gene in Pathway"] = 0  # Initialize, will be updated later
            to_be_update_pathway["gene_list"] = list()
            KOGatheringState = True
    del to_be_update_pathway, KOGatheringState
    pathways_df =  pd.DataFrame(pathways)
        
    return pathways_df


def _write_deg_foldchange_StrInfo(gene_list_df):
    """
    Writes the differential expression fold change information of genes in a string format.

    Parameters:
    gene_list_df (DataFrame): DataFrame containing gene expression information.

    Returns:
    str: String containing the gene ID, sub name, and fold change information of the top 5 genes.
    """
    gene_list_df = gene_list_df.sort_values(by="log2FoldChange", ascending=False)
    gene_list_df["fold_change"] = gene_list_df.log2FoldChange.apply(lambda x: 2 ** x)

    str_description = ""
    for idx, row in gene_list_df.reset_index().iterrows():
        if idx >= 5:
            break

        sub_name = ""

        if idx > 0:
            str_description += ", "
        str_description += f"{row.Geneid}({sub_name}{round(row['fold_change'], 2)})  "

    return str_description
    
def KEGG_Enrichment_Analysis(res_deg):
    """
    Perform KEGG enrichment analysis.

    Args:
        res_deg (pandas.DataFrame): DataFrame containing the differentially expressed genes (DEGs).

    Returns:
        pandas.DataFrame: DataFrame containing the results of the KEGG enrichment analysis.

    This function performs KEGG enrichment analysis by comparing the differentially expressed genes (DEGs) 
    with the gene lists associated with KEGG pathways. It calculates the number of DEGs and non-DEGs 
    that are present or absent in each pathway, and generates a summary table of the results.

    The function requires the `res_deg` DataFrame, which should contain the DEGs. It also relies on the 
    `_ConstructKEGGPathways` function to obtain the gene lists associated with KEGG pathways.

    Example usage:
        res = GetDEGs("OR-36h", "WT-36h", UpDown = "Down")
        kegg_result = KEGGAnalysis.KEGG_Enrichment_Analysis(res)
    """
    return_df = []
    pathways = _ConstructKEGGPathways()
    
    current_dir = os.path.dirname(os.path.abspath(__file__))
    total_gene_list = pd.read_csv(os.path.join(current_dir, "..", "Genomes", "Sclav_Gene_Info.csv"), index_col = 0).query("`gene_biotype` == 'protein_coding'").locus_tag.to_list()
    total_gene_list = set(total_gene_list)
    
    deg_list = set(res_deg.index.to_list())
    for idx, row in pathways.iterrows():
        pathway_gene_list = set(row.gene_list)
        
        n_DEG_Path = len(deg_list.intersection(pathway_gene_list))
        n_notDEG_Path = len(pathway_gene_list - deg_list)
        n_DEG_notPath = len(deg_list - pathway_gene_list)
        n_notDEG_notPath = len(total_gene_list - deg_list - pathway_gene_list)

        cur_pathway_table = pd.DataFrame([[n_notDEG_notPath, n_DEG_notPath],
                                          [n_notDEG_Path, n_DEG_Path]], columns = ['not DEG', 'DEG'], index = ['not in Pathway', 'in Pathway'])
        cur_pathway_table.loc['Total DEG'] = cur_pathway_table.sum(axis = 0)
        cur_pathway_table["Total pathway"] = cur_pathway_table.sum(axis = 1)

        assert cur_pathway_table.iloc[-1, -1] == len(total_gene_list)

        _, pval, _, _ = scipy.stats.chi2_contingency(cur_pathway_table.iloc[:-1, :-1])
        
        DEG_in_Path_info = _write_deg_foldchange_StrInfo(res_deg.query("Geneid in @deg_list.intersection(@pathway_gene_list)"))
        return_df.append( {"p value": pval,  "# of gene in Pathway": len(pathway_gene_list),  
                           "Unit": row["Unit"], "Subunit": row["Subunit"], "Pathway code": row["Pathway code"], "Pathway name": row["Pathway name"],
                           "# of DEGs in Pathways": n_DEG_Path,  "DEG in Pathways": DEG_in_Path_info})
        
    return pd.DataFrame(return_df).sort_values(by = "p value").reset_index(drop = True)