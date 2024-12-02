import pandas as pd

# Load the data from the files
taxonomy_df = pd.read_csv('taxonomy_fil_filtered.tsv', sep='\t')
feature_table_df = pd.read_csv('feature-table.tsv', sep='\t')

# Insert the Taxon_fil_filtered column as the second column in feature_table_df
feature_table_df.insert(1, 'Taxon_fil_filtered', taxonomy_df['Taxon_fil_filtered'])

# Save the updated feature table back to a file
feature_table_df.to_csv('feature-table_updated.tsv', sep='\t', index=False)

print("The 'Taxon_fil_filtered' column has been added as the second column in 'feature-table_updated.tsv'.")
