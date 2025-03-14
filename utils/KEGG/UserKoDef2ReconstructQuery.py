import os
import pandas as pd

user_ko_df = pd.read_csv("user_ko_definition.txt", sep = '\t', header = None).iloc[:, [0, 1]]
user_ko_df.columns = ["protein_id", "ko_number"]
user_ko_df.dropna(inplace = True)

sclav_df = pd.read_csv( os.path.join("..", "..", "Genomes", "Sclav_Gene_Info.csv")  )[["protein_id", "locus_tag"]]

result_df = []
for idx, row in user_ko_df.iterrows():
    cur_proteinid_df = sclav_df.query("protein_id == @row.protein_id")
    for locus_tag in cur_proteinid_df.locus_tag:
        result_df.append([locus_tag, row.ko_number])
result_df = pd.DataFrame(result_df, columns = ["locus_tag", "ko_number"])

result_df.to_csv("KEGGPathwayReconstructionQuery.txt", sep = '\t', index = False, header = False)