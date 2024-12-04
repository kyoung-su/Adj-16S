import pandas as pd
import os

# Set working directory
input_dir = '{set the directory}'

# Process all files ending with _matched_alignment_gaps.csv in the specified directory
for filename in os.listdir(input_dir):
    if filename.endswith('_matched_alignment_gaps.csv'):
        # Set file path
        file_path = os.path.join(input_dir, filename)
        
        # Read the CSV file
        df = pd.read_csv(file_path)
        
        # Extract only rows with even indices
        odd_rows_df = df.iloc[0::2]
        
        # Generate a new file name (e.g., original_matched_alignment_gaps.csv -> original_mag_extracted.csv)
        output_filename = filename.replace('_matched_alignment_gaps.csv', '_mag_extracted.csv')
        output_path = os.path.join(input_dir, output_filename)
        
        # Save the extracted data to a new CSV file
        odd_rows_df.to_csv(output_path, index=False)
        print(f"Processed and saved: {output_filename}")

print("All matching files have been processed and saved.")

