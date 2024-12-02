#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import os
import time

# Define the function to execute system commands
def run_system_commands():
    # Get the current working directory
    dir_file = os.getcwd()
    
    # List all files in the directory
    files = os.listdir(dir_file)
    
    # Find pairs of files ending with `_1_trim.fastp.fastq` and `_2_trim.fastp.fastq`
    file_pairs = []
    for file1 in files:
        if file1.endswith("_1_trim.fastp.fastq"):
            file2 = file1.replace("_1_trim.fastp.fastq", "_2_trim.fastp.fastq")
            if file2 in files:
                file_pairs.append((file1, file2))
    
    # Process each file pair
    for file1, file2 in file_pairs:
        # Construct the system command
        cmd = f"jtax.pl ref.fa 13primer.fa {file1} {file2} -join dj"
        print(f"Running: {cmd}")
        
        # Execute the command
        os.system(cmd)
        time.sleep(1)  # Pause to avoid overwhelming the system
            
# Run the function
if __name__ == "__main__":
    run_system_commands()
