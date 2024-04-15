#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Inputs: trrust_rawdata.human.tsv MAF(mutations.csv)

Outputs: tf - target (Activation or Repression)
        MAF of tfs with known regulation containing missense mutations
        some plots to explore
"""

import pandas as pd
import os
import sys
from Bio import SeqIO, Entrez
from io import StringIO
import requests
import matplotlib.pyplot as plt

pd.options.mode.chained_assignment = None  # default='warn'

def fetch_sequence_lengths(ids, id_type):
    sequence_lengths = {}
    
    # Choose the base URL based on the ID type
    base_url = "https://rest.uniprot.org/uniparc/" if id_type == 'uniparc' else "https://rest.uniprot.org/uniprotkb/"
    
    # ids have to be list
    for seq_id in ids:
        # Construct the URL for fetching the sequence data
        url = f"{base_url}{seq_id}.fasta"

        try:
            # Send a GET request to the UniProt or UniParc API
            response = requests.get(url)

            # Check if the request was successful
            if response.status_code == 200:
                # Parse the response content as FASTA data
                fasta_data = StringIO(response.text)
                records = SeqIO.parse(fasta_data, "fasta")

                # Calculate the length of the sequence from the FASTA content
                for record in records:
                    sequence_lengths[seq_id] = float(len(record.seq))
            else:
                print(f"Error fetching {seq_id}: Status Code {response.status_code}")
        except Exception as e:
            print(f"Error fetching {seq_id}: {str(e)}")
    
    return sequence_lengths

if len(sys.argv) != 3:
    print("Usage: program.py <trrust_rawdata.human.tsv> <mutations.csv> ")
else:
    filename = sys.argv[1]
    #base_filename = os.path.splitext(os.path.basename(filename))[0]
    maf = sys.argv[2]
    
    # Create a directory for saving results
    outputdir = "filtered_data"
    os.makedirs(outputdir, exist_ok=True)

############################     TFs filtering      ###########################


column_names = ['Gene1', 'Gene2', 'InteractionType', 'DOI']
df = pd.read_csv (filename, sep='\t', header=None, names=column_names )

# Remove tf pairs with unknown InteractionType
tf = df[df['InteractionType'] != 'Unknown']

# Group the DataFrame by 'Gene1' and 'Gene2'
grouped = tf.groupby(['Gene1', 'Gene2'])

# Filter out tf-target pairs with both 'Repression' and 'Activation' for the same target
for group_name, group_df in grouped:
    # Check if both 'Repression' and 'Activation' are in the 'InteractionType' values within the group
    if 'Repression' in group_df['InteractionType'].values and \
       'Activation' in group_df['InteractionType'].values:
        # Get the indices of rows to drop
        indices_to_drop = group_df.index
        tf = tf.drop(indices_to_drop)

# Reset the index of the filtered DataFrame
tf = tf.reset_index(drop=True)

# Save tf - target with known regulation
outfile1 = os.path.join(outputdir, 'trrust_processed.human.tsv')
tf.to_csv(outfile1, index=False)

print("TF - target with known interaction saved at", outputdir, 'trrust_processed.human.tsv \n')

tf_genes = tf["Gene1"].unique()

print("TFs were filtered./n")

############################     MAF Filtering      ###########################


MAF = pd.read_csv (maf, sep=',', low_memory=False)

# Filter MAF to keep only rows of the TFs selected and missense mutations.
MAF = MAF[MAF["Hugo_Symbol"].isin(tf_genes) & (MAF['Variant_Classification'] == 'Missense_Mutation') ]

# 1) Fetch the length for entries with no uniparc using swissprot id

# Extract the UniProt when the UniParc is NA and fetch its lengths
filtered_df = MAF[MAF['UNIPARC'].isna()]
uniprots = filtered_df['SWISSPROT'].tolist()
uniprots = list(set(uniprots))
lengths = fetch_sequence_lengths(uniprots , 'uniprot')
  
# Map UniProt IDs to sequence lengths using the lengths dictionary
filtered_df['Length'] = filtered_df['SWISSPROT'].map(lengths)

# 2) Fetch the length for entries that do have uniparc
uniparcs = MAF["UNIPARC"].tolist()
uniparcs = list(set(uniparcs))
lengths = fetch_sequence_lengths(uniparcs, 'uniparc')

# Map UNIPARC IDs to sequence lengths using the lengths dictionary
MAF['Length'] = MAF['UNIPARC'].map(lengths)

# Merge these 2 df by identifing the index values of the rows from filtered_df
filtered_indices = filtered_df.index.tolist()

# and update the original MAF DataFrame with the rows from filtered_df
MAF.loc[filtered_indices] = filtered_df

# count hoy many mutations one gene has
nr_mut_per_gene = MAF["Hugo_Symbol"].value_counts()

# Iterate through the rows of the MAF DataFrame
mutation_loads = list()
for index, row in MAF.iterrows():
    hugo_symbol = row['Hugo_Symbol']
    
    # Map Hugo symbol to mutation count
    mutation_count = nr_mut_per_gene.get(hugo_symbol)
    
    # Get the length value from the "Length" column
    length = row['Length']
    
    # Calculate the result and append it to the list
    result = mutation_count / int(length) if length != 0 else 0  # Avoid division by zero
    mutation_loads.append(result)

# Create a new column "Mutation_Count" in the MAF DataFrame
MAF['Mutation_load'] = mutation_loads

# Save the DataFrame as a CSV file in the output directory
df_filename = os.path.join(outputdir, 'TF_mutations_processed.csv')
MAF.to_csv(df_filename, index=False)



################          PLOTS         ###########################

# WRONG: Group the DataFrame by "Hugo_Symbol" and "UNIPARC" and take the max of "Mutation_Count" for each group
#grouped = MAF.groupby(["Hugo_Symbol", "UNIPARC"])["Mutation_load"].max().reset_index()

# CORRECTION : Group the DataFrame by "Hugo_Symbol" and take the max of "Mutation_load" for each group
grouped = MAF.groupby(["Hugo_Symbol"])["Mutation_load"].max().reset_index()


# Create a histogram of the mutational load
plt.hist(grouped['Mutation_load'], bins=100, edgecolor='k')
plt.xlabel('Mutation_load')
plt.ylabel('Frequency')
plt.title('Mutational Load Distribution')
plt.grid(True)

# Save the plot in the output directory
plot_filename = os.path.join(outputdir, 'mutational_load_histogram.png')
plt.savefig(plot_filename)

#######     Create a plot with the top 20 genes with highest mutational load      ###########

# Sort the DataFrame by mutational load in descending order
sorted_df = grouped.sort_values(by='Mutation_load', ascending=False)

# Select the top 20 rows from each group (UNIPARC and Hugo_Symbol combination)
top_20_groups = sorted_df.head(20).reset_index(drop=True)

# Create a bar plot for the top 20 rows
plt.figure(figsize=(10, 6))
plt.bar(
    range(len(top_20_groups)),
    top_20_groups['Mutation_load']
)
plt.xlabel('Hugo Name')
plt.ylabel('Mutational Load')
plt.title('Top 20 Genes with Highest Mutational Load')

# Set x-axis labels to Hugo_Symbol
plt.xticks(range(len(top_20_groups)), top_20_groups.apply(lambda row: f" ({row['Hugo_Symbol']})", axis=1), rotation=90)

plt.tight_layout()


# Save the plot in the output directory and the sorted df
plot_filename = os.path.join(outputdir, 'top_20_mutational_load_uniparc_hugo_symbol.png')
plt.savefig(plot_filename)


filename = os.path.join(outputdir, 'top_TFs.csv')
sorted_df.to_csv(filename, index=False)

plt.show()

print("Mutations were processed./n")
    

