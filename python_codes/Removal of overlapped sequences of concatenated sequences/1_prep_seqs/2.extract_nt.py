import os
from Bio import SeqIO

def extract_nucleotide_sequences():
    # Define the input and output directories
    current_directory = "{set the directory}"
    output_directory = "{set the directory}"
    
    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Loop through each file in the input directory
    for filename in os.listdir(current_directory):
        if filename.endswith("dj.fastq") or filename.endswith(".fasta"):
            input_path = os.path.join(current_directory, filename)
            output_path = os.path.join(output_directory, filename.replace("dj.fastq", "dj.fasta").replace(".fasta", "_nt.fasta"))

            with open(output_path, "w") as output_file:
                # Determine the file format based on the file extension
                file_format = "fastq" if filename.endswith(".fastq") else "fasta"
                
                # Read the file and save each sequence as .fasta
                for record in SeqIO.parse(input_path, file_format):
                    SeqIO.write(record, output_file, "fasta")

            print(f"Extracted sequences from {filename} -> {output_path}")

# Run the function
extract_nucleotide_sequences()

