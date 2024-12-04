from Bio import SeqIO
from Bio import AlignIO
import pandas as pd
import os
import subprocess
from glob import glob

# Define the output directory for saving CSV files
output_directory = "{set the directory}"

# Ensure the output directory exists
os.makedirs(output_directory, exist_ok=True)

# Function to extract base filename (up to .V13) for matching
def extract_base_name(filename):
    return filename.split(".V13")[0]

# Get list of FASTA files in the current directory
fasta_files = glob("*.fasta")

# Group files by base name
file_groups = {}
for filepath in fasta_files:
    base_name = extract_base_name(os.path.basename(filepath))
    if base_name not in file_groups:
        file_groups[base_name] = []
    file_groups[base_name].append(filepath)

# Process each matched file group
for base_name, files in file_groups.items():
    # Ensure we have exactly two files for each base name
    if len(files) == 2:
        # Identify dj and complement files based on naming
        dj_file = next((f for f in files if 'dj' in f), None)
        complement_file = next((f for f in files if 'complement' in f), None)

        if dj_file and complement_file:
            try:
                # Load sequences from FASTA files into dictionaries keyed by sequence IDs
                dj_sequences = {record.id: record for record in SeqIO.parse(dj_file, 'fasta')}
                complement_sequences = {record.id: record for record in SeqIO.parse(complement_file, 'fasta')}

                # Prepare to collect results
                results = []

                # Iterate over sequence IDs in the dj file and find matching IDs in the complement file
                for seq_id, dj_record in dj_sequences.items():
                    if seq_id in complement_sequences:
                        complement_record = complement_sequences[seq_id]

                        # Write each matched pair to a temporary FASTA file for alignment
                        with open('matched_pair.fasta', 'w') as match_file:
                            SeqIO.write([dj_record, complement_record], match_file, 'fasta')

                        # Perform pairwise alignment on the matched pair using MUSCLE
                        muscle_cmd = ["muscle", "-in", "matched_pair.fasta", "-out", "aligned_pair.fasta"]
                        subprocess.run(muscle_cmd, check=True)

                        # Read the aligned sequences
                        alignment = AlignIO.read('aligned_pair.fasta', 'fasta')

                        # Count gaps in each aligned sequence after the 300th nucleotide
                        for record in alignment:
                            gap_count = sum(1 for base in str(record.seq)[300:] if base == '-')
                            results.append({
                                'Sequence_ID': record.id,
                                'Gap_Count_After_300': gap_count,
                                'Aligned_Sequence': str(record.seq)
                            })

                # Create a DataFrame and save it to a CSV file in the output directory
                output_csv_name = os.path.join(output_directory, f"{base_name}_matched_alignment_gaps.csv")
                results_df = pd.DataFrame(results)
                results_df.to_csv(output_csv_name, index=False)

                print(f"Analysis of matched sequences completed for {base_name}. Results saved to '{output_csv_name}'.")

            except Exception as e:
                print(f"An error occurred while processing {base_name}: {e}")

            finally:
                # Clean up temporary files
                if os.path.exists('matched_pair.fasta'):
                    os.remove('matched_pair.fasta')
                if os.path.exists('aligned_pair.fasta'):
                    os.remove('aligned_pair.fasta')

