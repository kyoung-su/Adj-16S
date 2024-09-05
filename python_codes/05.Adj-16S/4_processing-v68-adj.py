import pandas as pd
import numpy as np


# Load the updated feature table
feature_table_df = pd.read_csv('feature-table_updated.tsv', sep='\t')

# Load the coefficient Excel file
coefficient_df = pd.read_excel('coefficient.xlsx')

# Create a mapping dictionary from the coefficient file
mapping_dict = pd.Series(coefficient_df.iloc[:, 1].values, index=coefficient_df.iloc[:, 0]).to_dict()

# Replace the names in the second column of feature_table_df using the mapping dictionary
feature_table_df.iloc[:, 1] = feature_table_df.iloc[:, 1].map(mapping_dict)


# Save the updated DataFrame to a new TSV file
feature_table_df.to_csv('feature-v68.tsv', sep='\t', index=False)

print("The names in the second column have been replaced and saved to 'feature-v68.tsv'.")

# Load the feature-v68.tsv file
df = pd.read_csv('feature-v68.tsv', sep='\t')

# Rename the second column's header
df.columns = [df.columns[0], 'coefficient-v68'] + list(df.columns[2:])

# Replace empty values and NaNs in the second column with 0.5
df['coefficient-v68'] = df['coefficient-v68'].replace('', np.nan)  # First, replace empty strings with NaN
df['coefficient-v68'] = df['coefficient-v68'].fillna(0.5)          # Then, replace NaNs with 0.5
df['coefficient-v68'] = df['coefficient-v68'].astype(float)        # Ensure the column is of type float

# Remove columns where any cell contains the substring 'V68'
df = df.loc[:, ~df.apply(lambda col: col.astype(str).str.contains('V13').any(), axis=0)]

# Shift the second column down by one row
df['coefficient-v68'] = df['coefficient-v68'].shift(1)

# Save the modified DataFrame to a new TSV file
df.to_csv('feature-v68.tsv', sep='\t', index=False)

print("The first row has been removed, the second column header has been renamed. The result is saved to 'feature-v68.tsv'.")
