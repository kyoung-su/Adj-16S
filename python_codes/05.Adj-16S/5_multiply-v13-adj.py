import pandas as pd

# Load the TSV file into a DataFrame
df = pd.read_csv('feature-v13.tsv', sep='\t')

# Make sure there are at least 3 columns in the DataFrame
if df.shape[1] < 3:
    raise ValueError("The input file must have at least 3 columns.")

# Store the first two rows (they will be skipped in calculations)
header_rows = df.iloc[:2].copy()

# Process the rest of the DataFrame
df_rest = df.iloc[2:].copy()

# Convert the columns of interest to numeric, coercing errors to NaN
df_rest.iloc[:, 1:] = df_rest.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')

# Perform the multiplication for the columns starting from the third one
def multiply_row(row):
    if pd.notna(row[1]):  # Check if the multiplier (second column value) is not NaN
        return row[2:] * row[1]
    else:
        return row[2:]

df_rest.iloc[:, 2:] = df_rest.apply(multiply_row, axis=1)

# Concatenate the header rows with the processed DataFrame
result_df = pd.concat([header_rows, df_rest], ignore_index=True)

# Save the modified DataFrame to a new TSV file
result_df.to_csv('feature-v13-adj.tsv', sep='\t', index=False)

print("The file has been processed and saved as 'feature-v13-adj.tsv'.")
