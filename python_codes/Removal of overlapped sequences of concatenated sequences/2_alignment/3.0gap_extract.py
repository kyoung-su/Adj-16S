import pandas as pd
import os

# Set working directory
input_dir = '{set the directory}'

# Process all files ending with _mag_extracted.csv in the specified directory
for filename in os.listdir(input_dir):
    if filename.endswith('_mag_extracted.csv'):
        # Set file path
        file_path = os.path.join(input_dir, filename)
        
        # Read the CSV file
        df = pd.read_csv(file_path)
        
        # Filter rows where 'Gap_Count_After_300' is 0
        filtered_df = df[df['Gap_Count_After_300'] == 0]
        
        # Generate a new text file name (e.g., original_mag_extracted.csv -> original_sequence_ids_with_zero_gaps.txt)
        output_filename = filename.replace('_mag_extracted.csv', '_sequence_ids_with_zero_gaps.txt')
        output_path = os.path.join(input_dir, output_filename)
        
        # Save the 'Sequence_ID' column to a text file
        filtered_df['Sequence_ID'].to_csv(output_path, index=False, header=False)
        print(f"Processed and saved: {output_filename}")

print("All matching files have been processed and saved.")

