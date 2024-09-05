#!/usr/bin/phtyon
# -*- coding: UTF-8 -*-

import os
import os.path
import time

class SYSTEM_CMD:

    def __init__(self):

        os.chdir("os.getcwd()")

        DIR_FILE = os.getcwd()
        
        for X in range(5000, 10000):
          
            FileName01 = f"SRR1300{X}_1_trim.fastp.fastq"
            FileName02 = f"SRR1300{X}_2_trim.fastp.fastq"
            
            if os.path.isfile(f"{DIR_FILE}{FileName01}") and os.path.isfile(f"{DIR_FILE}{FileName02}"):
                    CMD = f"jtax.pl ref.fa 13primer.fa {FileName01} {FileName02} -join dj"
                    os.system(CMD)
                    time.sleep(1)
                    
            else:
                pass

SYSTEM_CMD()
