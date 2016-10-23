# Montipora_Spawn_Timing


This repository includes data and analysis scripts to accompany:

Authors:  
Journal:   
Link:   
Description: 
This project examines gene expression during the normal spawning timing of corals that did and did not spawn. Samples were collected from 3 replicate _Montipora capitata_ colonies for each condition at 18:00, 20:00 and 00:00 the night of 6/6/16. RNA extracted from samples was prepped using Illumina TruSeq Sample Prep v2 Guide and sequenced 2x150bp on the Illumina MiSeq.
 
## Data Description

**Sample.ID** | **Date** | **Time** | **Spawn** | **Read.ID** | **File.ID** | **Read.Count**  
 ------ | ------ | ------ | ------ | ------ | ------ | ------ 
T4-1 | 6/6/16 | 18:00 | NO | Read1 | T4-1_S1_L001_R1_001.fastq.gz |
T4-1 | 6/6/16 | 18:00 | NO | Read2 | T4-1_S1_L001_R2_001.fastq.gz |
T4-6 | 6/6/16 | 18:00 | NO | Read1 | T4-6_S1_L001_R1_001.fastq.gz |
T4-6 | 6/6/16 | 18:00 | NO | Read2 | T4-6_S1_L001_R2_001.fastq.gz |
T4-8 | 6/6/16 | 18:00 | YES | Read1 | T4-8_S1_L001_R1_001.fastq.gz |
T4-8 | 6/6/16 | 18:00 | YES | Read2 | T4-8_S1_L001_R2_001.fastq.gz |
T4-10 | 6/6/16 | 18:00 | YES | Read1 | T4-10_S1_L001_R1_001.fastq.gz |
T4-10 | 6/6/16 | 18:00 | YES | Read2 | T4-10_S1_L001_R2_001.fastq.gz |
T4-16 | 6/6/16 | 18:00 | YES | Read1 | T4-16_S1_L001_R1_001.fastq.gz |
T4-16 | 6/6/16 | 18:00 | YES | Read2 | T4-16_S1_L001_R2_001.fastq.gz |
T4-17 | 6/6/16 | 18:00 | NO | Read1 | T4-17_S1_L001_R1_001.fastq.gz |
T4-17 | 6/6/16 | 18:00 | NO | Read2 | T4-17_S1_L001_R2_001.fastq.gz |
T5-1 | 6/6/16 | 20:00 | NO | Read1 | T5-1_S2_L001_R1_001.fastq.gz |
T5-1 | 6/6/16 | 20:00 | NO | Read2 | T5-1_S2_L001_R2_001.fastq.gz |
T5-6 | 6/6/16 | 20:00 | NO | Read1 | T5-6_S2_L001_R1_001.fastq.gz |
T5-6 | 6/6/16 | 20:00 | NO | Read2 | T5-6_S2_L001_R2_001.fastq.gz |
T5-8 | 6/6/16 | 20:00 | YES | Read1 | T5-8_S2_L001_R1_001.fastq.gz |
T5-8 | 6/6/16 | 20:00 | YES | Read2 | T5-8_S2_L001_R2_001.fastq.gz |
T5-10 | 6/6/16 | 20:00 | YES | Read1 | T5-10_S2_L001_R1_001.fastq.gz |
T5-10 | 6/6/16 | 20:00 | YES | Read2 | T5-10_S2_L001_R2_001.fastq.gz |
T5-16 | 6/6/16 | 20:00 | YES | Read1 | T5-16_S2_L001_R1_001.fastq.gz |
T5-16 | 6/6/16 | 20:00 | YES | Read2 | T5-16_S2_L001_R2_001.fastq.gz |
T5-17 | 6/6/16 | 20:00 | NO | Read1 | T5-17_S2_L001_R1_001.fastq.gz |
T5-17 | 6/6/16 | 20:00 | NO | Read2 | T5-17_S2_L001_R2_001.fastq.gz |
T7-1 | 6/7/16 | 00:00 | NO | Read1 | T7-1_S3_L001_R1_001.fastq.gz |
T7-1 | 6/7/16 | 00:00 | NO | Read2 | T7-1_S3_L001_R2_001.fastq.gz |
T7-6 | 6/7/16 | 00:00 | NO | Read1 | T7-6_S3_L001_R1_001.fastq.gz |
T7-6 | 6/7/16 | 00:00 | NO | Read2 | T7-6_S3_L001_R2_001.fastq.gz |
T7-8 | 6/7/16 | 00:00 | YES | Read1 | T7-8_S3_L001_R1_001.fastq.gz |
T7-8 | 6/7/16 | 00:00 | YES | Read2 | T7-8_S3_L001_R2_001.fastq.gz |
T7-10 | 6/7/16 | 00:00 | YES | Read1 | T7-10_S3_L001_R1_001.fastq.gz |
T7-10 | 6/7/16 | 00:00 | YES | Read2 | T7-10_S3_L001_R2_001.fastq.gz |
T7-16 | 6/7/16 | 00:00 | YES | Read1 | T7-16_S3_L001_R1_001.fastq.gz |
T7-16 | 6/7/16 | 00:00 | YES | Read2 | T7-16_S3_L001_R2_001.fastq.gz |
T7-17 | 6/7/16 | 00:00 | NO | Read1 | T7-17_S3_L001_R1_001.fastq.gz |
T7-17 | 6/7/16 | 00:00 | NO | Read2 | T7-17_S3_L001_R2_001.fastq.gz |



6 lanes with 3 samples per lane for library prep and sequencing 2 x 150 on MiSeq at HIMB
Samples submitted: 20160629
Estimated Seq Oct 2016
see [**Library_Prep_Info.csv**]()
see [**Sample.Info.csv**]()

## Notes about Files

# Directory Structure


## notebooks:
Jupyter Notebook files for bioinformatic analysis

## Processing:
Quality and quantity results from RNA extraction and library preparation 

## Protocols:
Sample processing protocols for RNA extraction and library preparation

## RAnalysis:
	### Data
	### Output
	### Scripts

## SEQ_Data:
	### Clean_Data
	### Cleaned_QC_Files
	### Raw_Data
	### Raw_QC_Files
 











