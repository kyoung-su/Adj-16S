import os
from Bio import SeqIO, Seq
import pandas as pd

# Define paths for input and output directories
zero_gap_dir = '{set the directory}'
fastq_dir = '{set the directory}'
output_dir = '{set the directory}'

# Find pairs of matching files
zero_gap_files = [f for f in os.listdir(zero_gap_dir) if f.endswith('_sequence_ids_with_zero_gaps.txt')]
fastq_files = [f for f in os.listdir(fastq_dir) if f.endswith('.V13_1_trim.fastp.dj.fastq')]

# Process each pair of matching files
for zero_gap_file in zero_gap_files:
    # Extract the base name to find the matching FASTQ file
    base_name = zero_gap_file.replace('_sequence_ids_with_zero_gaps.txt', '')
    fastq_file = f"{base_name}.V13_1_trim.fastp.dj.fastq"
    
    if fastq_file in fastq_files:
        # Load the list of Sequence IDs to analyze
        with open(os.path.join(zero_gap_dir, zero_gap_file), 'r') as f:
            target_ids = set(line.strip() for line in f if line.strip())
        
        # Load sequences from the FASTQ file and filter to include only target IDs
        sequences = [record for record in SeqIO.parse(os.path.join(fastq_dir, fastq_file), 'fastq') if record.id in target_ids]
        
        # Prepare a list to store results
        results = []
        
        # Iterate over each target sequence in the FASTQ file
        for record in sequences:
            sequence = str(record.seq)
            sequence_id = record.id
            quality = record.letter_annotations["phred_quality"]  # Quality scores associated with the sequence

            # Track positions of matched sequences
            matched_ranges = []

            # Iterate through range expansions for the comparison
            for length in range(1, 207):  # Expanding from 1 bp up to 207 bp
                start1 = 300 - length   # Starting position for first range (backwards expansion)
                end1 = 300              # Ending position for first range

                start2 = 301               # Starting position for second range
                end2 = 301 + length        # Ending position for second range (forwards expansion)

                # Ensure indices are within sequence bounds
                if start1 >= 0 and end2 <= len(sequence):
                    # Extract sequences for current range
                    range1_seq = sequence[start1:end1]
                    range2_seq = sequence[start2:end2]

                    # Check if the sequences match
                    if range1_seq == range2_seq:
                        matched_ranges.append((start1, end1, start2, end2))  # Store matching ranges

            # Remove nucleotides and quality scores corresponding to one of the matching ranges
            if matched_ranges:
                # Use the first matched range to remove nucleotides
                start1, end1, start2, end2 = matched_ranges[0]
                
                # Create new sequence and quality with one instance of the duplicate removed
                new_sequence = (sequence[:start1] + 
                                sequence[start1:end1] +  # Keep the nucleotides between the first and second range
                                sequence[end2:])  # Keep the rest after the second range
                
                new_quality = (quality[:start1] +
                               quality[start1:end1] +  # Keep the quality scores for kept nucleotides
                               quality[end2:])  # Keep the quality scores after the second range
            else:
                # No matches found, keep the original sequence and quality
                new_sequence = sequence
                new_quality = quality

            # Append results with the modified sequence and quality
            results.append({
                'Sequence_ID': sequence_id,
                'Final_Sequence': new_sequence,
                'Final_Quality': new_quality
            })

        # Convert results to a DataFrame and save to CSV
        output_csv = os.path.join(output_dir, f"{base_name}_nucleotide_matching_results_filtered.csv")
        results_df = pd.DataFrame(results)
        results_df.to_csv(output_csv, index=False)

        # Write modified sequences to a new FASTQ file
        output_fastq = os.path.join(output_dir, f"{base_name}.V13_dj_trim.fastq")
        with open(output_fastq, 'w') as output_fastq_file:
            for result in results:
                seq_record = SeqIO.SeqRecord(
                    Seq.Seq(result['Final_Sequence']),
                    id=result['Sequence_ID'],
                    description=''
                )
                seq_record.letter_annotations["phred_quality"] = result['Final_Quality']
                SeqIO.write(seq_record, output_fastq_file, 'fastq')

        print(f"Comparison completed for {fastq_file}. Results saved to '{output_csv}' and trimmed FASTQ to '{output_fastq}'.")

