import pandas as pd

# Reading the TSV File:
df = pd.read_csv('taxonomy.tsv', sep='\t')

# Filtering the Taxon Column:
df['Taxon_filtered'] = df['Taxon'].str.split('; g__').str[0]

# Saving the First Filtered Data:
df.to_csv('taxonomy_filtered.tsv', sep='\t', index=False)

print("Filtered taxonomy data saved as 'taxonomy_filtered.tsv'")

import pandas as pd

# Reading the Filtered TSV File:
df = pd.read_csv('taxonomy_filtered.tsv', sep='\t')

# Further Filtering the Taxon_filtered Column:
df['Taxon_fil_filtered'] = df['Taxon_filtered'].str.split('; f__').str[1]

# Saving the Further Filtered Data:
df.to_csv('taxonomy_fil_filtered.tsv', sep='\t', index=False)

print("Filtered taxonomy data saved as 'taxonomy_fil_filtered.tsv'")
