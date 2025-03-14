import pandas as pd
import scipy.stats
import os
import re

current_dir = os.path.dirname(os.path.abspath(__file__))

def _LoadModuleDF():
    with open( os.path.join(current_dir, "KEGG", "ModuleReconstructionResult.txt"), "r") as f:
        lines = f.readlines()
    lines.append("M00000 END MODULE (0)\n")
        
    module_df = []
    module_code = 'empty'
    for idx in range(len(lines)):
        line = lines[idx].strip()
        
        pattern = r'^(M\d{5}|K\d{5}|CRV15_RS\d+)'
        if not bool(re.match(pattern, line)):   # if its text
            if not bool(re.match(pattern, lines[idx+1].strip())):
                large_category = line
            else:
                small_category = line
        
        elif bool(re.match(r'^(M\d{5})', line)):  # if its module
            # check if `module_code` is defined
            if module_code != 'empty':
                module_df.append({"large_category": large_category, "small_category": small_category, 
                                  "module_code": module_code, "module_name": module_name, "N_components": int(N_components), "N_genes": len(gene_list), "gene_list": gene_list})
            
            match = re.search(r'^(M\d{5})\s(.*?)(?=\s+\(\d+\))\s\((\d+)\)', line)   
            if match:
                module_code, module_name, N_components = match.groups()
                gene_list = []
            else:
                raise ValueError("Unrecognized Module line")
            
        elif bool(re.match(r'^(CRV15_RS\d+)', line)):
            gene_list.extend(line.split(", "))
        
        elif bool(re.match(r'^(K\d{5})', line)):
            continue
        
        else:
            raise ValueError("Unrecognized line")
        
    module_df = pd.DataFrame(module_df)
    module_df = module_df[(module_df['N_components'] >= 3) & (module_df['N_genes'] >= 4)].reset_index(drop=True)
    return module_df


def _LoadKEGGPathwayDF():

    total_gene_list = pd.read_csv(os.path.join(current_dir, "..", "Genomes", "Sclav_Gene_Info.csv"), index_col = 0).query("`gene_biotype` == 'protein_coding'").locus_tag.to_list()
    total_gene_list = set(total_gene_list)    
    
    with open( os.path.join(current_dir, "KEGG", "PathwayReconstruction_Curated.txt"), "r") as f:
        lines = f.readlines()
    lines.append("00000 END MODULE (0)\n")
        
    pathway_df = []
    pathway_code = 'empty'
    for idx in range(len(lines)):
        
        line = lines[idx].strip()
                
        pattern = r'^(\d{5}|K\d{5}|CRV15_RS\d+)'
        if not bool(re.match(pattern, line)):   # if its text
            if not bool(re.match(pattern, lines[idx+1].strip())):
                large_category = line
            else:
                small_category = line
        
        elif bool(re.match(r'^(\d{5})', line)):  # if its pathway
            # check if `pathway_code` is defined
            if pathway_code != 'empty':
                gene_list = list(set(gene_list).intersection(total_gene_list))
                pathway_df.append({"large_category": large_category, "small_category": small_category, 
                                   "Pathway_code": pathway_code, "Pathway_name": pathway_name, 
                                   "N_components": int(N_components), "N_genes": len(gene_list), "gene_list": gene_list})
            
            match = re.search(r'^(\d{5})\s(.*?)(?=\s+\(\d+\))\s\((\d+)\)', line)   
            if match:
                pathway_code, pathway_name, _ = match.groups()
                N_components = 0
                gene_list = []
            else:
                print(line)
                raise ValueError("Unrecognized Module line")
            
        elif bool(re.match(r'^(CRV15_RS\d+)', line)):
            N_components += 1
            gene_list.extend(line.split(", "))
        
        elif bool(re.match(r'^(K\d{5})', line)):
            continue
        
        else:
            raise ValueError("Unrecognized line")
        
    pathway_df = pd.DataFrame(pathway_df)
    pathway_df = pathway_df[(pathway_df['N_components'] >= 3) & (pathway_df['N_genes'] >= 4)].reset_index(drop=True)
    return pathway_df


def ConstructOntologyDF():
    module_df = _LoadModuleDF()
    pathway_df = _LoadKEGGPathwayDF()
    
    # Rename the columns in module_df and pathway_df so they match
    module_df = module_df.rename(columns={"module_code": "KEGG code", "module_name": "KEGGterm name"})
    pathway_df = pathway_df.rename(columns={"Pathway_code": "KEGG code", "Pathway_name": "KEGGterm name"})

    # Concatenate the two DataFrames
    ontology_df = pd.concat([module_df, pathway_df], ignore_index=True)
    return ontology_df


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
    pathways = ConstructOntologyDF()
    
    total_gene_list = pd.read_csv(os.path.join(current_dir, "..", "Genomes", "Sclav_Gene_Info.csv"), index_col = 0).query("`gene_biotype` == 'protein_coding'").locus_tag.to_list()
    total_gene_list = set(total_gene_list)
    
    deg_list = set(res_deg.index.to_list()) if isinstance(res_deg, pd.DataFrame) else set(res_deg)
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
        
        DEG_in_Path_info = _write_deg_foldchange_StrInfo(res_deg.query("Geneid in @deg_list.intersection(@pathway_gene_list)"))   \
                                if isinstance(res_deg, pd.DataFrame)    \
                                else ", ".join(list(deg_list.intersection(pathway_gene_list)))
        return_df.append( {"p value": pval,  "# of gene in Pathway": len(pathway_gene_list), 
                           "Large category": row["large_category"], "Small category": row["small_category"],
                           "KEGG code": row["KEGG code"], "KEGGterm name": row["KEGGterm name"],
                           "# of DEGs in Pathways": n_DEG_Path,  "DEG in Pathways": DEG_in_Path_info})
        
    return pd.DataFrame(return_df).sort_values(by = "p value").reset_index(drop = True)