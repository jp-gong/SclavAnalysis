import os
import pandas as pd
from tqdm import tqdm
from Bio import SeqIO

def gff_df_file(file_path, leave_only_gene_info=True):
    """
    Reads a GFF file and returns a DataFrame with gene-related rows.
    """
    gff_df = pd.read_csv(
        file_path,
        comment="#",
        sep="\t",
        names=["contig", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    )
    if leave_only_gene_info:
        for feature in ["region", "direct_repeat", "riboswitch", "exon"]:
            gff_df = gff_df[gff_df["feature"] != feature].reset_index(drop=True)
        gff_df.drop(columns=["score", "frame"], inplace=True)
    return gff_df

def Get_geneattr_from_gff(gff_df: pd.DataFrame):
    """
    Parses the attribute fields of paired GFF rows to extract gene attributes.
    """
    def parse_attribute(attr1, attr2):
        attr1_dict = {att.split("=")[0] + "_1": att.split("=")[1] for att in attr1.split(";")}
        attr2_dict = {att.split("=")[0] + "_2": att.split("=")[1] for att in attr2.split(";")}
        attr1_dict.update(attr2_dict)
        return attr1_dict

    gene_attr_list = []
    # Process rows in pairs
    for i in tqdm(range(0, len(gff_df), 2), desc="Processing GFF rows"):
        row1 = gff_df.loc[i].copy()
        row2 = gff_df.loc[i + 1].copy()
        # Ensure that the paired rows have matching start and end coordinates
        assert (row1.start == row2.start) and (row1.end == row2.end)
        attr = parse_attribute(row1.attribute, row2.attribute)
        row1.drop(["source", "feature", "attribute"], inplace=True)
        for key, value in attr.items():
            row1[key] = value
        row1["contig"] = "chromosome" if row1["contig"] == "NZ_CP027858.1" else "plasmid"
        if "gene_synonym_1" in row1:
            row1["gene_1"] += " " + row1["gene_synonym_1"]
            row1.drop("gene_synonym_1", inplace=True)
        if "Ontology_term_2" in row1:
            total_string = ""
            for go_category in ["go_function_2", "go_process_2", "go_component_2"]:
                if go_category in row1:
                    for go_term in row1[go_category].split(","):
                        total_string += f"; {go_term.split('|')[0]} [GO:{go_term.split('|')[1]}]"
            row1["gene_ontology"] = total_string[2:].replace("%2C", ",")
        gene_attr_list.append(row1)
    
    gene_attr_df = pd.DataFrame(gene_attr_list).reset_index(drop=True)

    # Drop unnecessary columns
    cols_to_drop = [
        "gbkey_1", "ID_1", "transl_table_2", "gene_2", "pseudo_2", "partial_2",
        "anticodon_2", "Name_1", "Parent_2", "locus_tag_2", "gbkey_2",
        "Name_2", "ID_2", "Dbxref_2", "old_locus_tag_1", "Note_2",
        "partial_1", "pseudo_1", "start_range_1", "start_range_2", "end_range_1", "end_range_2",
        "Ontology_term_2", "go_function_2", "go_process_2", "go_component_2"
    ]
    gene_attr_df.drop(columns=cols_to_drop, inplace=True, errors='ignore')
    
    gene_attr_df.rename(columns={
        "gene_biotype_1": "gene_biotype",
        "locus_tag_1": "locus_tag",
        "inference_2": "inference",
        "product_2": "product",
        "protein_id_2": "protein_id",
        "gene_1": "gene_name"
    }, inplace=True)
    
    # Rearrange columns
    gene_attr_df = gene_attr_df[[
        "contig", "locus_tag", "gene_name", "start", "end", "strand",
        "gene_biotype", "product", "protein_id", "gene_ontology"
    ]]
    return gene_attr_df

def clean_uniprot_mapping(file_path):
    """
    Reads and cleans the Uniprot mapping file,
    returning a DataFrame filtered for Streptomyces clavuligerus.
    """
    uniprot_df = pd.read_csv(file_path, sep="\t")
    uniprot_df = uniprot_df.query("Organism == 'Streptomyces clavuligerus'").reset_index(drop=True)
    
    uniprot_cleaned_list = []
    for protein_id, group in uniprot_df.groupby('From'):
        cur_gene_list = []
        for gene in group["Gene Names"]:
            if pd.isna(gene):
                continue
            cur_gene_list.extend(gene.split())
        cur_row = {
            "protein_id": protein_id,
            "Uniprot Entry": " ".join(group["Entry"]),
            "Uniprot Protein names": group["Protein names"].iloc[0],
            "Uniprot Gene Names": " ".join(sorted(set(cur_gene_list)))
        }
        uniprot_cleaned_list.append(cur_row)
    return pd.DataFrame(uniprot_cleaned_list)

def main():
    # Process the GFF file to generate the gene attribute DataFrame
    gff_file = os.path.join("reference_genome", "genomic.gff")
    gff_df = gff_df_file(gff_file)
    gene_attr_df = Get_geneattr_from_gff(gff_df)
    
    # Process the Uniprot mapping file
    uniprot_file = os.path.join("data", "uniprot_refsefAccession_mapping.tsv")
    uniprot_cleaned_df = clean_uniprot_mapping(uniprot_file)
    
    # Merge gene attributes with Uniprot mapping on protein_id
    gene_annot_df = gene_attr_df.merge(uniprot_cleaned_df, how='outer', on='protein_id')
    assert len(gene_annot_df) == len(gene_attr_df)
    
    # Save the merged DataFrame to CSV
    output_file = "Sclav_Gene_Info.csv"
    gene_annot_df.to_csv(output_file, index=False)
    print(f"Generated {output_file}")

if __name__ == "__main__":
    main()
