import os
from Bio import SeqIO

def reverse_complement_fastq_to_fasta():
    # Define the directory where the .fastq files are located
    current_directory = "{set the directory}"
    
    # Loop through each file in the directory
    for filename in os.listdir(current_directory):
        if filename.endswith(".fastq"):
            input_path = os.path.join(current_directory, filename)
            output_path = os.path.join(current_directory, filename.replace(".fastq", "_reverse_complement.fasta"))

            with open(output_path, "w") as output_file:
                # Read .fastq file and create reverse complement
                for record in SeqIO.parse(input_path, "fastq"):
                    # Generate the reverse complement of the sequence
                    reverse_complement_seq = record.seq.reverse_complement()
                    # Update the record with the reverse complement sequence and description
                    record.seq = reverse_complement_seq
                    record.description = record.id + " reverse complement"
                    # Write to output as .fasta format
                    SeqIO.write(record, output_file, "fasta")

            print(f"Processed {filename} -> {output_path}")

# Run the function
reverse_complement_fastq_to_fasta()

