import subprocess
import os

# List of scripts to execute
scripts = [
    "1_extract_taxa.py",
    "2_dataframe_rearrange.py",
    "3_processing-v13-adj.py",
    "4_processing-v68-adj.py",
    "5_multiply-v13-adj.py",
    "6_multiply-v68-adj.py",
    "7_sum-adj.py"
]

# Create a result directory if it does not exist
result_dir = "result"
os.makedirs(result_dir, exist_ok=True)

# Execute each script in sequence
for script in scripts:
    print(f"Running {script}...")
    
    # Run the script and capture output
    result = subprocess.run(["python", script], capture_output=True, text=True)
    
    # Define output file path
    output_file = os.path.join(result_dir, f"{os.path.splitext(script)[0]}.log")
    
    # Write the output to the result directory
    with open(output_file, "w") as file:
        file.write(result.stdout)
    
    # Write the error to the result directory if any occurs
    if result.returncode != 0:
        error_file = os.path.join(result_dir, f"{os.path.splitext(script)[0]}_error.log")
        with open(error_file, "w") as file:
            file.write(result.stderr)
        print(f"Error occurred while running {script}, see {error_file} for details.")
        break  # Stop further execution if an error occurs

    print(f"Output written to {output_file}")
