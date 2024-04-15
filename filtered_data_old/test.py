#!/usr/bin/env python3

# Some UNIPARCs are NAs so we have to group only by HUGO otherwise genes like TP53(na UNIPARC) will not 
# be included in the grouped df 

import matplotlib
matplotlib.use('Agg')  # Set the backend before importing pyplot
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("TF_mutations_processed.csv", sep=',', low_memory=False)
# Group the DataFrame by "Hugo_Symbol" and "UNIPARC" and take the max of "Mutation_Count" for each group
grouped = df.groupby(["Hugo_Symbol"])["Mutation_load"].max().reset_index()

# Create a histogram of the mutational load
plt.hist(grouped['Mutation_load'], bins=100, edgecolor='k')
plt.xlabel('Mutation_load')
plt.ylabel('Frequency')
plt.title('Mutational Load Distribution')
plt.grid(True)

# Save the plot in the output directory

plt.savefig('test.png')




