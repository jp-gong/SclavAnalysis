import pandas as pd
import scipy.stats
import os

def Build_goterms_from_gbff():
    """
    Builds Gene Ontology (GO) terms from a gbff file.

    Args:
        total_gene_list (list): List of genes to consider.

    Returns:
        pandas.DataFrame: DataFrame containing the GO terms with their associated information.
    """
    current_dir = os.path.dirname(os.path.abspath(__file__))
    total_gene_list = pd.read_csv(os.path.join(current_dir, "..", "Genomes", "Sclav_Gene_Info.csv"), index_col = 0).query("`gene_biotype` == 'protein_coding'").locus_tag.to_list()
    total_gene_list = set(total_gene_list)
        
    # Read `gbff` file
    with open(os.path.join(current_dir, "../Genomes/reference_genome/genomic.gbff")) as f:
        text = f.readlines()
    chrm, plsmd = text[76:127303], text[239860:271607]
    text = chrm + plsmd

    # Remove spaces at the beginning of each line
    for line in text:
        assert line[:5] == "     "
    text = [line[5:] for line in text]

    # Make `go_text` DataFrame
    locus_tag_flag = 0
    go_text = []
    for idx, line in enumerate(text):
        if line[0] != ' ':
            keyword = line[:16].strip()
            locus_tag_flag = 0
            continue
        
        line = line[16:-1]
        
        if line.startswith("/locus_tag"):
            locus_tag = line[12:-1]
            locus_tag_flag = 1
        if line.startswith("/GO"):
            while text[idx+1][16] != '/':
                line += ' ' + text[idx+1][16:-1]
                idx += 1
            assert locus_tag_flag == 1
            go_text.append({"gene": locus_tag, "text": line})
    go_text = pd.DataFrame(go_text)

    # parse `go_text` to `go_terms`
    go_terms = dict()
    for idx, (gene, line) in go_text.iterrows():
        assert len(line.split('"')) == 3
        domain = line.split('"')[0][4:-1]
        
        line = line.split('"')[1]
        assert len(line.split(' - ')) == 2        
        go_term, explanation = line.split(' - ')
        
        if gene not in total_gene_list:
            continue
        
        if go_term not in go_terms:
            go_terms[go_term] = {"domain": domain, "go_term": go_term, "name": explanation, "gene_list": [gene]}
        else:
            go_terms[go_term]["gene_list"].append(gene)
            assert go_terms[go_term]["name"] == explanation
            
    go_terms = pd.DataFrame(go_terms).T
    go_terms.insert(3, "# of gene in GO", go_terms["gene_list"].apply(len))
    go_terms = go_terms.sort_values(by = "# of gene in GO", ascending = False).reset_index(drop = True)
    return go_terms



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



def GO_Enrichment_Analysis(res_deg):
    """
    Perform Gene Ontology (GO) enrichment analysis.

    Parameters:
    res_deg (DataFrame): DataFrame containing the differentially expressed genes (DEGs).

    Returns:
    DataFrame: DataFrame containing the results of the GO enrichment analysis, sorted by p-value.

    """
    
    return_df = []
    go_terms = Build_goterms_from_gbff()
    
    current_dir = os.path.dirname(os.path.abspath(__file__))
    total_gene_list = pd.read_csv(os.path.join(current_dir, "..", "Genomes", "Sclav_Gene_Info.csv"), index_col = 0).query("`gene_biotype` == 'protein_coding'").locus_tag.to_list()
    total_gene_list = set(total_gene_list)
    
    deg_list = set(res_deg.index.to_list())
    for idx, row in go_terms.iterrows():
        goterm_gene_list = set(row.gene_list)
        
        n_DEG_Path = len(deg_list.intersection(goterm_gene_list))
        n_notDEG_Path = len(goterm_gene_list - deg_list)
        n_DEG_notPath = len(deg_list - goterm_gene_list)
        n_notDEG_notPath = len(total_gene_list - deg_list - goterm_gene_list)

        cur_pathway_table = pd.DataFrame([[n_notDEG_notPath, n_DEG_notPath],
                                          [n_notDEG_Path, n_DEG_Path]], columns = ['not DEG', 'DEG'], index = ['not in Pathway', 'in Pathway'])
        cur_pathway_table.loc['Total DEG'] = cur_pathway_table.sum(axis = 0)
        cur_pathway_table["Total pathway"] = cur_pathway_table.sum(axis = 1)

        assert cur_pathway_table.iloc[-1, -1] == len(total_gene_list)

        # display(cur_pathway_table)

        _, pval, _, _ = scipy.stats.chi2_contingency(cur_pathway_table.iloc[:-1, :-1])
        
        DEG_in_Path_info = _write_deg_foldchange_StrInfo(res_deg.query("Geneid in @deg_list.intersection(@goterm_gene_list)"))
        return_df.append( {"p value": pval,  "# of gene in GO": len(goterm_gene_list),  
                           "domain": row["domain"], "go_term": row["go_term"],  "GO name": row["name"].replace(" [Evidence IEA]", ""),
                           "# of DEGs in GO": n_DEG_Path,  "DEG in pathway": DEG_in_Path_info})
        
    return pd.DataFrame(return_df).sort_values(by = "p value").reset_index(drop = True)