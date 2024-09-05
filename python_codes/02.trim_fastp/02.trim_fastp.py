#!/usr/bin/phtyon
# -*- coding: UTF-8 -*-

import os
import os.path
import time

class SYSTEM_CMD:

    def __init__(self):

        os.chdir(os.getcwd())

        DIR_FILE = os.getcwd()

        MAX_NUM = 9000
        
        for X in range(5000, MAX_NUM):

         
        	FileName01 = f"SRR1300{X}_1_fastp.fastq"
        	FileName02 = f"SRR1300{X}_2_fastp.fastq"
        	FileName03 = f"SRR1300{X}_1_trim.fastp.fastq"
        	FileName04 = f"SRR1300{X}_2_trim.fastp.fastq"
        	
        	if os.path.isfile(f"{DIR_FILE}{FileName01}") and os.path.isfile(f"{DIR_FILE}{FileName02}"):
        	    CMD = f"fastp -i {FileName01} -I {FileName02} -b 300 -B 208 -o {FileName03} -O {FileName04}"
        	    os.system(CMD)
        	    time.sleep(1)
        	    
        	else:
        	    pass

SYSTEM_CMD()
