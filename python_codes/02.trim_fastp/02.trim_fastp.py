#!/usr/bin/python
# -*- coding: UTF-8 -*-

import os
import time
import glob

class SYSTEM_CMD:
    def __init__(self):
        os.chdir(os.getcwd())
        DIR_FILE = os.getcwd()

        # Find all files ending with _1.fastp.fastq and _2.fastp.fastq
        file_pairs = {}
        for file in glob.glob("*_1.fastp.fastq"):
            base_name = file.replace("_1.fastp.fastq", "")
            paired_file = f"{base_name}_2.fastp.fastq"
            if os.path.isfile(paired_file):
                file_pairs[base_name] = (file, paired_file)

        for base_name, (FileName01, FileName02) in file_pairs.items():
            FileName03 = f"{base_name}_1_trim.fastp.fastq"
            FileName04 = f"{base_name}_2_trim.fastp.fastq"
            
            CMD = f"fastp -i {FileName01} -I {FileName02} -b 300 -B 208 -o {FileName03} -O {FileName04}"
            os.system(CMD)
            time.sleep(1)

SYSTEM_CMD()
#If the trimming position changes, the user only needs to modify the number.
