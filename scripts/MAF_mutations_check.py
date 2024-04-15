import sys
import os
import pandas as pd
import numpy as np


if len(sys.argv) != 4:
    print("Usage: program.py <TF_mutations.csv> <metatable.csv> <GENE> ")
else:
    maf = sys.argv[1]
    cancermuts = sys.argv[2]
    GENE = sys.argv[3]
    
    # Create a directory for saving results
    outputdir = "../mavisp_missing_mutations"
    os.makedirs(outputdir, exist_ok=True)

output_columns = ['name', 'site', 'type', 'function', 'reference']
output_df = pd.DataFrame(columns=output_columns)


maf = pd.read_csv (maf, sep=',', low_memory=False)
maf = maf[maf['Hugo_Symbol'] == GENE ]
maf = maf.dropna(subset=['HGVSp', 'HGVSp_Short'])


cancermuts = pd.read_csv(cancermuts, sep= ',')

#maf_mutations = maf['HGVSp', 'HGVSp_short']
maf_mutations = maf[["HGVSp" , "HGVSp_Short"]].astype(str) # p.Gln165AlafsTer6. , pY65T

lst = list()
cancermuts = cancermuts.dropna(subset=["alt_aa"])


for index, row in cancermuts.iterrows():
    cancermut = str(row["ref_aa"]) + str(row["aa_position"])  +  str(row["alt_aa"] )
    lst.append(cancermut)

    

for index, row in maf_mutations.iterrows():
    if row['HGVSp_Short'][2:] not in lst:

        # Append mutation information as a new row to the output DataFrame
        output_df.loc[len(output_df)] = {'name': row['HGVSp_Short'][2:],
                                      'site': row['HGVSp'],
                                      'type': 'mutation',
                                      'function': '',
                                      'reference': ''}


# Save to CSV or use the DataFrame as needed
output_df.to_csv(outputdir +'/missing_mutuations.csv', sep=';', index=False)

print("MAF file contained " ,len(maf_mutations), "mutations.")
print(len(output_df), "mutation(s) aren't covered from cancermuts dataset.Import them manually. File saved at missing_mutuations.csv")




