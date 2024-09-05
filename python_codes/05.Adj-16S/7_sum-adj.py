import pandas as pd
import numpy as np

# Load the TSV files into DataFrames
df_v13 = pd.read_csv('feature-v13-adj.tsv', sep='\t')
df_v68 = pd.read_csv('feature-v68-adj.tsv', sep='\t')

# Check if the number of rows and columns match
if df_v13.shape != df_v68.shape:
    raise ValueError("The dimensions of the two files do not match.")

# Create a new DataFrame for the result
result_df = df_v13.copy()

# Define a function to safely add and round numeric values
def safe_add_round(val1, val2):
    try:
        # Convert values to numeric, forcing errors to NaN
        num1 = pd.to_numeric(val1, errors='coerce')
        num2 = pd.to_numeric(val2, errors='coerce')
        # Add and round if both are numeric
        if not pd.isna(num1) and not pd.isna(num2):
            return round(num1 + num2)
        else:
            return np.nan  # or some other placeholder for missing/invalid data
    except:
        return np.nan

# Perform element-wise addition and rounding for all columns
for i in range(1, len(df_v13.columns)):
    result_df.iloc[:, i] = [
        safe_add_round(val1, val2) for val1, val2 in zip(df_v13.iloc[:, i], df_v68.iloc[:, i])
    ]

# Remove the second column (index 1)
result_df = result_df.drop(result_df.columns[1], axis=1)

# Save the result to a new TSV file
result_df.to_csv('feature-adj.tsv', sep='\t', index=False)

print("Feature adjustment completed and saved to 'feature-adj.tsv'.")
