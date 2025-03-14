import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import joblib

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
RefGenome_DIR = os.path.join(BASE_DIR, "..", "Genomes", "reference_genome")
strains = ["WT", "C1", "OR", "NL", "MA"]

def _load_translate_table():
    translate_table = dict()
    start_codons = dict()
    stop_codons = dict()
    
    cds_seq_dict = {record.id: str(record.seq) for record in SeqIO.parse(os.path.join(RefGenome_DIR, "cds_from_genomic.fna"), "fasta")}
    protein_dict = {record.id: str(record.seq) for record in SeqIO.parse(os.path.join(RefGenome_DIR, "protein.faa"), "fasta")}

    for i, cds_seq in cds_seq_dict.items():
        j = i.split("_cds_")[1].rsplit("_", 1)[0]
        
        if not j.startswith("WP"):
            continue
        
        assert ( (len(cds_seq) / 3 - 1) == len(protein_dict[j]) ), "ERROR"
        
        # Iterate through the CDS sequence in steps of 3 (codon size) and the protein sequence
        for k in range(0, len(cds_seq)-3, 3):  # -3 to adjust for the stop codon not being translated
            codon = cds_seq[k:k+3]
            amino_acid = protein_dict[j][k//3]  # Convert k from nucleotide to protein position
            
            if k == 0:   # Start codon is always translated to Methionine
                assert amino_acid == "M", "ERROR"
                if codon not in start_codons:
                    start_codons[codon] = 0
                start_codons[codon] += 1
                continue
            
            # Check if codon is already in translate_table with a different amino acid
            if codon in translate_table and translate_table[codon] != amino_acid:
                print(f"Conflict detected for codon {codon}: {translate_table[codon]} vs {amino_acid},  range = {k}")
                raise ValueError(f"Conflict detected for codon {codon}: {translate_table[codon]} vs {amino_acid}")
            
            # Assign the codon to the amino acid in the translation table
            translate_table[codon] = amino_acid
            
        # Check the stop codon
        stop_codon = cds_seq[-3:]
        if stop_codon not in stop_codons:
            stop_codons[stop_codon] = 0
        stop_codons[stop_codon] += 1
    
    return translate_table, start_codons, stop_codons

def _load_mutation_data(strains):
    vcf_data = {}
    for strain in strains:
        vcf_file_path = os.path.join(BASE_DIR, "..", "MutationResults", f"{strain}_variants.vcf")
        vcf_data[strain] = pd.read_csv(vcf_file_path, sep="\t", comment="#",
                                       names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "bam"])
        vcf_data[strain] = vcf_data[strain].drop(columns=["ID", "FILTER", "FORMAT", "bam"])
        vcf_data[strain]["CHROM"] = vcf_data[strain]["CHROM"].apply(lambda x: "Linear" if x == "NZ_CP027858.1" else "Plasmid")
    return vcf_data 

def _get_amino_acid_change(gene_df, gene_df_row, mut_row, protein_dict):
    if gene_df_row['rel.POS'].startswith("-"):
        return "promoter region"
    
    if gene_df_row['protein_id'] in protein_dict:
        if gene_df_row['rel.POS'] in ["+1", "+2", "+3"]:
            return "start codon mut"
        
        ref_len = len(mut_row["REF"])
        alt_len = len(mut_row["ALT"])
        if ref_len != alt_len:
            if (ref_len - alt_len) % 3 != 0:
                return "frameshift"
            else:
                return "inframe"
        
        protein_seq = protein_dict[gene_df_row['protein_id']]
        return "CDS mutation"
    else:
        return ""

def _search_gene_by_pos(gene_df, contig, row):
    pos = row["POS"]
    protein_dict = {record.id: str(record.seq) for record in SeqIO.parse(os.path.join(RefGenome_DIR, "protein.faa"), "fasta")}
    
    def calculate_rel_pos(row, pos):
        if row["strand"] == "+":    return f'+{pos - row["start"] + 1}' if pos - row["start"] >= 0 else f'-{row["start"] - pos}'
        else:                       return f'+{row["end"] - pos + 1}' if row["end"] - pos >= 0 else f'-{pos - row["end"]}'

    if contig.lower() in ["linear", "NZ_CP027858.1"]:     contig = "chromosome"
    if contig.lower() in ["plasmid", "NZ_CP027859.1"]:    contig = "plasmid"
    
    cur_df = gene_df[gene_df["contig"] == contig.lower()]
    assert len(cur_df) != 0, f"Contig {contig} not found"

    cur_df = cur_df[(cur_df["start"] - 40 <= pos) & (cur_df["end"] >= pos) & (cur_df["strand"] == "+") |
                    (cur_df["start"] <= pos) & (cur_df["end"] + 40 >= pos) & (cur_df["strand"] == "-")]
    
    do_not_search = len(cur_df) == 0
    if do_not_search:
        cur_df.loc[0] = ""
    
    cur_df["POS"] = pos
    if not do_not_search:
        cur_df["rel.POS"] = cur_df.apply(lambda x: calculate_rel_pos(x, pos), axis=1)
    cur_df["REF"] = row["REF"]
    cur_df["ALT"] = row["ALT"]
    cur_df["QUAL"] = row["QUAL"]
    cur_df["contig"] = row["CHROM"]
    
    if not do_not_search:
        cur_df["Mutation"] = cur_df.apply(lambda x: _get_amino_acid_change(gene_df, x, row, protein_dict), axis=1)
    
    return cur_df


def _get_genes_from_mut_info(gene_df, strain):
    cur_vcf_df = _load_mutation_data(strains)[strain]
    
    dfs = []
    for idx, row in cur_vcf_df.iterrows():
        temp_df = _search_gene_by_pos(gene_df, row["CHROM"], row)
        if not temp_df.empty:
            temp_df["OverlappedGene"] = len(temp_df) > 1
            dfs.append(temp_df)
        
    dfs = pd.concat(dfs, ignore_index=True)
    dfs = dfs.dropna(subset=["Mutation"]).reset_index(drop=True)
            
    return dfs


def LoadMutGeneInfo():
    os.makedirs(os.path.join(BASE_DIR, "MutationCacheData"), exist_ok=True)
    MutGeneInfo_path = os.path.join(BASE_DIR, "MutationCacheData", "MutGeneInfo.joblib")
    if os.path.exists(MutGeneInfo_path):
        MutGeneInfo = joblib.load(MutGeneInfo_path)
    else:
        SclavGeneInfo_path = os.path.join(BASE_DIR, "..", "Genomes", "Sclav_Gene_Info.csv")
        gene_df = pd.read_csv(SclavGeneInfo_path)
        MutGeneInfo = {strain: _get_genes_from_mut_info(gene_df, strain) for strain in strains}
        joblib.dump(MutGeneInfo, MutGeneInfo_path)
    return MutGeneInfo

def GetMutatedGenesFromGeneList(gene_list):
    MutGeneInfo = LoadMutGeneInfo()
    
    # Dictionary to hold mutated genes for each strain filtered by gene list
    mutated_genes = {}
    for strain in strains:
        mutated_genes[strain] = MutGeneInfo[strain][MutGeneInfo[strain]["locus_tag"].isin(gene_list)].copy()
        
    # List to store dataframes for each strain after adding a new column for strain        
    all_data_frames = []
    for strain, df in mutated_genes.items():
        df['Strain'] = strain  # Add a new column for strain
        all_data_frames.append(df)
    combined_df = pd.concat(all_data_frames, ignore_index=True).drop(columns = ['QUAL'])
    
    # Group the DataFrame by all columns except 'Strain' to find duplicates and aggregate the 'Strain' column,
    # joining the values with ', '. The 'as_index=False' makes the grouped columns regular columns in the resulting dataframe.
    combined_df = combined_df.groupby(by = combined_df.columns.difference(['Strain']).tolist(), as_index = False, dropna = False).agg({'Strain': ', '.join})
        
    return combined_df

def determine_mutation_type(nucleotide_mutation_df):
    nucleotide_complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    translate_table, start_codons, stop_codons = _load_translate_table()
    cds_seq_dict = {record.id.split("_cds_")[1].rsplit("_", 1)[0]: str(record.seq) for record in SeqIO.parse(os.path.join(RefGenome_DIR, "cds_from_genomic.fna"), "fasta")}
    protein_dict = {record.id: str(record.seq) for record in SeqIO.parse(os.path.join(RefGenome_DIR, "protein.faa"), "fasta")}
    
    working_df = nucleotide_mutation_df.copy()
    working_df['rel.POS'] = working_df['rel.POS'].astype(int)
    working_df['rel.AminoAcidPOS'] = np.where((working_df['gene_biotype'] == 'protein_coding') & (working_df['Mutation'] == 'CDS mutation'),
                                              (working_df['rel.POS'] - 1) // 3 + 1, working_df['rel.POS'])
    
    protein_mut_df = []
    for (cur_locus_tag, AminoAcidPos, cur_MutationType), cur_df in working_df.groupby(['locus_tag', 'rel.AminoAcidPOS', 'Mutation']):
        assert cur_df['gene_biotype'].nunique() == 1, "Error: Multiple gene biotypes for the same gene."
        
        # if is not a protein coding gene
        cur_gene_biotype = cur_df['gene_biotype'].iloc[0]
        if cur_gene_biotype != 'protein_coding':
            assert len(cur_df) == 1, f"Error: Multiple mutations for the same gene."
            protein_mut_df.append({'contig': cur_df['contig'].iloc[0], 'locus_tag': cur_locus_tag, 'gene_name': cur_df['gene_name'].iloc[0], 'gene_biotype': cur_gene_biotype,
                                   'protein_id': cur_protein_id, 'product': cur_df['product'].iloc[0],
                                   'REF N.T.': cur_df['REF'].iloc[0], 'ALT N.T.': cur_df['ALT'].iloc[0], 'rel.POS': AminoAcidPos,
                                   'AminoAcidPos': '-', 'REF A.A.': '-', 'ALT A.A.': '-', 'Mutation Type': cur_MutationType})
            continue
        
        # if is a protein coding gene
        assert cur_df['protein_id'].nunique() == 1, "Error: Multiple genes with the same protein_id."
        assert cur_df['strand'].nunique() == 1, "Error: Multiple strands for the same gene."
        cur_protein_id = cur_df['protein_id'].iloc[0]
        cur_strand = cur_df['strand'].iloc[0]
        
        assert cur_df['Mutation'].nunique() == 1, "Error: Multiple types of mutations at the same position."
        cur_MutationType = cur_df['Mutation'].iloc[0]
        if cur_MutationType == 'CDS mutation':
            AminoAcidPos = int(AminoAcidPos)
            original_codon = cds_seq_dict[cur_protein_id][(AminoAcidPos-1)*3:AminoAcidPos*3]

            altered_codon = [original_codon[0], original_codon[1], original_codon[2]]
            for _, row in cur_df.sort_values(by = 'rel.POS').iterrows():
                rel_POS_in_codon = (row['rel.POS']-1) % 3
                cur_REF = row['REF'] if cur_strand == '+' else nucleotide_complement[row['REF']]
                cur_ALT = row['ALT'] if cur_strand == '+' else nucleotide_complement[row['ALT']]
                assert altered_codon[rel_POS_in_codon] == cur_REF, "Error: Reference nucleotide does not match the sequence."
                altered_codon[rel_POS_in_codon] = cur_ALT
            altered_codon = ''.join(altered_codon)
            
            # if original codon is stop codon
            if (AminoAcidPos-1) == len(protein_dict[cur_protein_id]):
                if (original_codon not in translate_table) and (altered_codon not in translate_table):
                    NewMutationType = 'Nonsense'
                    cur_AminoAcid = '*'
                    altered_AminoAcid = '*'
                elif (original_codon not in translate_table) and (altered_codon in translate_table):
                    NewMutationType = 'NonStop'
                    cur_AminoAcid = '*'
                    altered_AminoAcid = translate_table[altered_codon]
                else:
                    raise ValueError("Error: Stop codon is not translated correctly.")
                
            # if original codon is not stop codon
            else:            
                cur_AminoAcid = protein_dict[cur_protein_id][AminoAcidPos-1]
                assert translate_table[original_codon] == cur_AminoAcid, "Error: Codon does not match the amino acid."
                
                if altered_codon not in translate_table:
                    altered_AminoAcid = '*'
                    NewMutationType = 'Nonsense'
                else:
                    altered_AminoAcid = translate_table[altered_codon]
                    NewMutationType = 'Silent' if altered_AminoAcid == cur_AminoAcid else 'Missense'
            
            protein_mut_df.append({'contig': cur_df['contig'].iloc[0], 'locus_tag': cur_locus_tag, 'gene_name': cur_df['gene_name'].iloc[0], 'gene_biotype': cur_gene_biotype,
                                   'protein_id': cur_protein_id, 'product': cur_df['product'].iloc[0],
                                   'REF N.T.': original_codon, 'ALT N.T.': altered_codon, 'rel.POS': '-',
                                   'AminoAcidPos': AminoAcidPos, 'REF A.A.': cur_AminoAcid, 'ALT A.A.': altered_AminoAcid, 'Mutation Type': NewMutationType})
        
        # if is not a CDS mutation (ex. promoter region, frameshift, ...)
        else:
            assert len(cur_df) == 1, f"Error: Multiple mutations for the same gene."
            protein_mut_df.append({'contig': cur_df['contig'].iloc[0], 'locus_tag': cur_locus_tag, 'gene_name': cur_df['gene_name'].iloc[0], 'gene_biotype': cur_gene_biotype,
                                   'protein_id': cur_protein_id, 'product': cur_df['product'].iloc[0],
                                   'REF N.T.': cur_df['REF'].iloc[0], 'ALT N.T.': cur_df['ALT'].iloc[0], 'rel.POS': AminoAcidPos,
                                   'AminoAcidPos': '-', 'REF A.A.': '-', 'ALT A.A.': '-', 'Mutation Type': cur_MutationType})
        
    protein_mut_df = pd.DataFrame(protein_mut_df)
    return protein_mut_df

def GetMutatedGenes():
    MutGeneInfo = LoadMutGeneInfo()
    
    def DetailsDescription(cur_df):
        context = ""
        
        for MutType, cur_type_df in cur_df.groupby("Mutation Type"):
            context += f"{MutType}("
            for idx, row in cur_type_df.iterrows():
                if MutType in ["Missense", "Nonsense", "NonStop", "Silent"]:
                    context += f"{row['REF A.A.']}{row['AminoAcidPos']}{row['ALT A.A.']},"
                else:
                    context += f"{row['rel.POS']}{row['REF N.T.'].lower()}>{row['ALT N.T.'].lower()},"
            context = context[:-1] + "), "
            
        context = context[:-2]
        return context
    
    # Dictionary to hold mutated genes for each strain
    mutated_genes = {}
    for strain in strains:
        mutated_genes[strain] = determine_mutation_type(MutGeneInfo[strain])
            
    for strain in strains:
        CleanedMutatedGenesInfo = []
        for locus_tag, cur_df in mutated_genes[strain].groupby('locus_tag'):
            MutTypes = set(cur_df["Mutation Type"].unique())
            # remove "Silent" element from MutTypes set
            if "Silent" in MutTypes:
                MutTypes.remove("Silent")
            
            if len(MutTypes) == 0:
                to_append_row = {'contig': cur_df['contig'].iloc[0], "locus_tag": locus_tag, "gene_name": cur_df["gene_name"].iloc[0], "gene_biotype": cur_df["gene_biotype"].iloc[0],
                                 "protein_id": cur_df["protein_id"].iloc[0], "product": cur_df["product"].iloc[0],
                                 "Mutation Type": "Silent", "Details": DetailsDescription(cur_df)}
            
            elif len(MutTypes) == 1:
                to_append_row = {'contig': cur_df['contig'].iloc[0], "locus_tag": locus_tag, "gene_name": cur_df["gene_name"].iloc[0], "gene_biotype": cur_df["gene_biotype"].iloc[0],
                                 "protein_id": cur_df["protein_id"].iloc[0], "product": cur_df["product"].iloc[0],
                                 "Mutation Type": MutTypes.pop(), "Details": DetailsDescription(cur_df)}
            else:
                to_append_row = {'contig': cur_df['contig'].iloc[0], "locus_tag": locus_tag, "gene_name": cur_df["gene_name"].iloc[0], "gene_biotype": cur_df["gene_biotype"].iloc[0],
                                 "protein_id": cur_df["protein_id"].iloc[0], "product": cur_df["product"].iloc[0],
                                 "Mutation Type": "Multiple", "Details": DetailsDescription(cur_df)}
            CleanedMutatedGenesInfo.append(to_append_row)
        mutated_genes[strain] = pd.DataFrame(CleanedMutatedGenesInfo)
    
    return mutated_genes  # Dictionary of DataFrames : {strain: mutated_genes_df}

def GetCommonMutatedGenes(mutated_genes: dict, strain1: str, strain2: str):
    strain1_mutated_genes = mutated_genes[strain1]
    strain2_mutated_genes = mutated_genes[strain2]
    
    assert strain1_mutated_genes.columns.equals(strain2_mutated_genes.columns), "Columns of the two DataFrames are not the same."
    
    result_df = pd.merge(strain1_mutated_genes, strain2_mutated_genes, on=list(strain1_mutated_genes.columns), how='inner', suffixes=('_' + strain1, '_' + strain2))
    return result_df

def GetDiffMutatedGenes(mutated_genes: dict, strain1: str, strain2: str):  # Get df1 - df2
    df1 = mutated_genes[strain1]
    df2 = mutated_genes[strain2]
    
    assert df1.columns.equals(df2.columns), "Columns of the two DataFrames are not the same."
    
    compare_columns = [col for col in df1.columns if col not in ['Mutation Type', 'Details']]
    result_df = pd.merge(df1, df2, on=compare_columns, how='outer', indicator=True, suffixes=('_df1', '_df2'))
    result_df = result_df[result_df['_merge'] == 'left_only'].drop(columns=['_merge'])
    
    result_df = result_df[ (result_df['Mutation Type_df1'] != result_df['Mutation Type_df2']) | (result_df['Details_df1'] != result_df['Details_df2']) ]
    result_df.rename(columns={'Mutation Type_df1': f'{strain1} Mutation Type', 'Details_df1': f'{strain1} Details',
                              'Mutation Type_df2': f'{strain2} Mutation Type', 'Details_df2': f'{strain2} Details'}, inplace=True)
    
    return result_df
    

if __name__ == "__main__":
    LoadMutGeneInfo()