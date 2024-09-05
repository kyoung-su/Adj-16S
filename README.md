# Adj-16S
Enhancement of Gut Microbial Diversity and Functional Profiling

**Features of Adj-16S**

1. Confirm Phred quality score (Q score)
To truncate sequence length, we first confirm the position at which the median Q score drops below 20 for both forward and reverse sequences. 

2. Concatenating paired-end sequences 
After we obtained length trimmed paired-end sequences, we concatenated forward and reverse sequences.

3. Use correction coefficient values
Utilizing mock datasets (Zymo, ZIEL-I, ZIEL-II) and the SILVA database, we calculated the correction coefficient values for each family using weighted averages, covering 22 families. Weighted coefficient values were calculated by using measured relative abundance obtained from V1-V3 and V6-V8 regions. Users can understand more easily after you read graphical representation of “workflows of Adj-16S”. 
(you can download mock community dataset in https://doi.org/10.1128/msphere.01202-20) (SRP291583)

4. Enhancing microbiome diversity and functional profiling
Raw data (phyloseq form) of adjusted relative abundance using Adj-16S method can be used for diverse R packages and PICRUSt2, not only mock community dataset, but you can also apply this methodology to diverse datasets which have paired-end sequences. Our 16S rRNA amplicon datasets (derived from human gut sample) using MiSeq platform are available on Sequence Read Archive (SRA), with project ID: PRJNA1088906.

**Workflows of Adj-16S**

![image](https://github.com/user-attachments/assets/3b8b3358-c8d4-4e73-9705-878ca92bdc5a)


Our developing Adj-16S method based on V1-V3 and V6-V8 regions is explained in detail in the “Advancing gut microbiota profiling accuracy with correction coefficient-based adjustments for dual 16S rRNA reads” section of the paper.
We will guide majorly how to use Adj-16S method in this document. Deatiled description about the calculation of correction coefficient values will be confirimed in our coming paper (You can confirm equations about the correction coefficient values for each family.)
Workflow figure was partially created using BioRender.com.

**Requirements**

**Software**

fastp (https://github.com/OpenGene/fastp) (version >= 0.23.2)

QIIME2 (https://github.com/qiime2/qiime2) (version >= 2022.2)

JTax (https://github.com/TLlab/JTax)

Python (version >= 3.7)

R (https://www.r-project.org/)

pandas (https://pandas.pydata.org/pandas-docs/version/0.11.0/install.html)

phyloseq (https://github.com/joey711/phyloseq) (optional)

PICRUSt2 (https://github.com/picrust/picrust2) (optional)


**16S rRNA raw data**

Users can download 16S rRNA raw data (SRP291583 and PRJNA1088906) from NCBI (https://www.ncbi.nlm.nih.gov/)


**16S rRNA Databases** 

SILVA (https://www.arb-silva.de/)

Greengenes2 (https://github.com/biocore/q2-greengenes2)

RDP (https://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/)

Users can download and align the sequences to diverse 16S rRNA databases.


**Other**

Operating system (Linux or Mac)


**Usage**

#1. Trimming sequences using fastp
We provide customized python code (**01.fastp.py**). In the script, you can easily change the code considering your file name. (Installation and guide about fastp, follow the above github link)

#2. Check Q score of paired-end sequences using QIIME2
Confirm the Q score of paired-end reads with QIIME2 platform. Before merging or concatenating paired-end sequences, we need to confirm where the median Q score falls below Q20.

#3. Make concatenated sequences using JTax 
To obtain merged sequences, you just follow QIIME2 document as described. But, QIIME2 platform cannot concatenate sequences. So, we used JTax platform to concatenate forward and reverse sequences. 
We can obtain concatenated sequences with direct-joining and inside-out methods. We also provide the customized python codes (**02.trim_fastp.py** and **03. JTax.py**).
First, you can trim sequences with appropriate positions which can be fixed in #2 using QIIME2. We provide customized python code 02, to trim easily for user’s availability. You can adjust (Line XX of code 02) trimming positions and file names.
Second, you can concatenate sequences with two methods (DJ and IO) using python code 03. After then, you can obtain the concatenated sequence files. In this step, you have to perform these python codes to V1-V3 and V6-V8 regions, separately. (As amplified sequence lengths are different, we cannot apply codes to these sequences at once). 
Third, perform a python **codes 04** (Clear_string_N.py and Clear_string_I.py). After concatenating sequences, ‘NNNNNNNN’ and ‘IIIIIIII’ are placed in the connected position between forward and reverse reads. These gaps can interfere with the next step of analysis and have to be removed. 

#4. Apply the correction coefficient values for each microbiome sample
First, users can easily merge V13-derived and V68-derived table.qza files with “qiime feature-table merge” command line (Please find a detailed command line with QIIME2 document). Second, users can obtain feature-table.tsv from table.qza using “biom convert -i feature-table.biom -o feature-table.tsv --to-tsv”. Third, perform taxonomic classification of V13-DJ and V68-DJ derived ASVs by searching 16S rRNA database such as SILVA, RDP, and greengenes2 using command line in QIIME2 document. Users can obtain taxonomy.tsv from taxonomy.qza.

1_extract_taxa.py and 2_dataframe_rearrange.py: Load family names from taxonomy.tsv and place them to the out-id column in feature-table.tsv. 
3_processing-v13-adj.py and 4_processing-v68-adj.py: Change each correction coefficient value to the corresponding family using coefficient.xlsx file. 
5_multiply-v13-adj.py and 6_multiply-v68-adj.py: Multiply the ASV value by the correction coefficient value and save it as feature-v13-adj.tsv and feature-v68-adj.tsv.
7_sum-adj.py: Apply this code to sum ASVs of feature-v13-adj.tsv and feature-v68-adj.tsv. Finally, users can obtain feature-adj.tsv file. Change only the columns of the feature-adj.tsv file to the same name as the sample name to be used in the metadata, save it, and then use the feature-adj.tsv file with phyloseq or PICRUSt2.

Users can perform these codes separately or run them all at once with the **05.Adj-16S.py** code.

#5 Use diverse R packages (phyloseq, microbiomemarker)
Using the feature-adj.tsv file, users can obtain more accurate microbiome composition and PICRUSt2-based functional profile data. Additionally, users can load the ASVs data with phyloseq (R package) and analyze microbiome data such as alpha and beta diversity.


**Applications**

To obtain taxonomic and functional profiling and perform statistical analyses, you can apply produced sequences for diverse R packages and PICRUSt2.
R packages (phyloseq, microbiomemarker, reshape2, tidyverse)
PICRUSt2 (https://github.com/picrust/picrust2)


**References**

1. Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 34, i884-i890 (2018).
2. Bolyen E, et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol 37, 852-857 (2019).
3. Liu T, et al. Joining Illumina paired-end reads for classifying phylogenetic marker sequences. BMC Bioinformatics 21, 105 (2020).
4. McMurdie PJ, Holmes S. phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. PLoS One 8, e61217 (2013).
5. Douglas GM, et al. PICRUSt2 for prediction of metagenome functions. Nat Biotechnol 38, 685-688 (2020).
