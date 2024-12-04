import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Define paths for input and output directories
trimmed_dir = '{set the directory}'
original_dir = '{set the directory}'
output_dir = '{set the directory}'

# Find pairs of matching files
trimmed_files = [f for f in os.listdir(trimmed_dir) if f.endswith('_dj_trim.fastq')]
original_files = [f for f in os.listdir(original_dir) if f.endswith('_1_trim.fastp.dj.fastq')]

# Process each pair of matching files
for trimmed_file in trimmed_files:
    # Extract the base name to find the matching original file
    base_name = trimmed_file.replace('_dj_trim.fastq', '')
    original_file = f"{base_name}_1_trim.fastp.dj.fastq"
    
    if original_file in original_files:
        # Load the trimmed sequences into a dictionary
        trimmed_sequences = {
            record.id: record for record in SeqIO.parse(os.path.join(trimmed_dir, trimmed_file), "fastq")
        }
        
        # Define output file path
        output_file = os.path.join(output_dir, f"{base_name}_1_dj_trim_final.fastq")
        
        # Open the original file and replace sequences if the ID matches
        with open(output_file, "w") as output_handle:
            for record in SeqIO.parse(os.path.join(original_dir, original_file), "fastq"):
                if record.id in trimmed_sequences:
                    # Create a new SeqRecord with trimmed sequence and quality
                    trimmed_record = trimmed_sequences[record.id]
                    new_record = SeqRecord(
                        trimmed_record.seq,
                        id=record.id,
                        description=record.description,
                        letter_annotations={"phred_quality": trimmed_record.letter_annotations["phred_quality"]}
                    )
                    SeqIO.write(new_record, output_handle, "fastq")
                else:
                    # Write the original record if no match is found
                    SeqIO.write(record, output_handle, "fastq")
        
        print(f"Replacement completed for {original_file}. Saved as {output_file}.")

