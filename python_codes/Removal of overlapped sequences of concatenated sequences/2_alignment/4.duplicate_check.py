from Bio import SeqIO
import pandas as pd
import os

# Set input and output directories
input_dir = '{set the directory}'
output_dir = os.path.join(input_dir, 'nucleotide_matching_results')

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Process all files ending with _1_trim.fastp.dj_nt.fasta in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith('_1_trim.fastp.dj_nt.fasta'):
        # Set file path
        file_path = os.path.join(input_dir, filename)
        
        # Load sequence file
        sequences = list(SeqIO.parse(file_path, 'fasta'))
        results = []

        # Perform analysis on each sequence
        for record in sequences:
            sequence = str(record.seq)
            sequence_id = record.id
            match_found = False

            # Expand comparison range (1 bp ~ 208 bp)
            for length in range(1, 207):
                start1 = 300 - length  # Start position of the first range
                end1 = 300             # End position of the first range

                start2 = 301           # Start position of the second range
                end2 = 301 + length    # End position of the second range

                # Check if the range is within sequence bounds
                if start1 >= 0 and end2 <= len(sequence):
                    range1_seq = sequence[start1:end1]
                    range2_seq = sequence[start2:end2]

                    # Check for matching sequences
                    if range1_seq == range2_seq:
                        results.append({
                            'Sequence_ID': sequence_id,
                            'Range1_Start': start1 + 1,  # Use 1-based index
                            'Range1_End': end1,
                            'Range2_Start': start2,
                            'Range2_End': end2,
                            'Matching_Length': length
                        })
                        match_found = True
                        break  # Stop expanding if a matching sequence is found

            # If no matching sequence is found
            if not match_found:
                results.append({
                    'Sequence_ID': sequence_id,
                    'Range1_Start': None,
                    'Range1_End': None,
                    'Range2_Start': None,
                    'Range2_End': None,
                    'Matching_Length': 0
                })

        # Convert results to DataFrame and save as CSV
        output_filename = filename.replace('_1_trim.fastp.dj_nt.fasta', '_nucleotide_matching_results.csv')
        output_path = os.path.join(output_dir, output_filename)
        results_df = pd.DataFrame(results)
        results_df.to_csv(output_path, index=False)
        print(f"Processed and saved: {output_filename}")

print("All matching files have been processed and saved.")

