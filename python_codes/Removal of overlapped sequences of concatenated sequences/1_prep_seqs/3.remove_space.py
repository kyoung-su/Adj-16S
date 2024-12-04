import os

# Function to process fasta file
def process_fasta_file(input_path, output_path):
    with open(input_path, 'r') as file:
        lines = file.readlines()

    # Process each line to remove text after the first space on lines that start with '>'
    processed_lines = []
    for line in lines:
        if line.startswith('>'):
            line = line.split()[0] + '\n'  # Keep only the part before the first space
        processed_lines.append(line)

    # Write the modified lines to the output file
    with open(output_path, 'w') as file:
        file.writelines(processed_lines)

# Define input and output directories
input_directory = "{set the directory}"
output_directory = "{set the directory}"

# Ensure the output directory exists
os.makedirs(output_directory, exist_ok=True)

# Process each '.nt.fasta' file in the input directory
for filename in os.listdir(input_directory):
    if filename.endswith('_nt.fasta'):
        input_path = os.path.join(input_directory, filename)
        output_path = os.path.join(output_directory, filename)
        
        # Process and save the modified file
        process_fasta_file(input_path, output_path)

print("Processing complete for all '.nt.fasta' files.")

