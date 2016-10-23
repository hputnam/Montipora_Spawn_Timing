# Analysis for Montipora_Spawn_Timing
Data uploaded and analyzed on Poire 
ssh hputnam@galaxy.geodata.hawaii.edu

# Upload data to server
```scp -r /Users/hputnam/MyProjects/Montipora_Spawn_Timing/Data/H_Putnam_RNA_run1-32840854 hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Data```

# Used the following adapter files to make barcodes file
>TruSeq_universal_F
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
> Genomic_DNA_oligonucleotide_sequences_Adapters_F
GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG
> Genomic_DNA_oligonucleotide_sequences_Adapters_R
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
> Genomic_DNA_Sequencing_Primer
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>TruSeq_idx_Adapter_AR001
GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR002
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR003
GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR004
GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR005
GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR006
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR007
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR008
GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR009
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR010
GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR011
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR012
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR013
GATCGGAAGAGCACACGTCTGAACTCCAGTCAAGTCAAAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR014
GATCGGAAGAGCACACGTCTGAACTCCAGTCAAGTTCCAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR015
GATCGGAAGAGCACACGTCTGAACTCCAGTCAATGTCAAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR016
GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTCCAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR018
GATCGGAAGAGCACACGTCTGAACTCCAGTCAGTCCGCAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR019
GATCGGAAGAGCACACGTCTGAACTCCAGTCAGTGAAAAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR020
GATCGGAAGAGCACACGTCTGAACTCCAGTCAGTGGCCAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR021
GATCGGAAGAGCACACGTCTGAACTCCAGTCAGTTTCGAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR022
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTACGAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR023
GATCGGAAGAGCACACGTCTGAACTCCAGTCAGAGTGGAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR025
GATCGGAAGAGCACACGTCTGAACTCCAGTCAACTGATAATCTCGTATGCCGTCTTCTGCTTG
>TruSeq_idx_Adapter_AR027
GATCGGAAGAGCACACGTCTGAACTCCAGTCAATTCCTAATCTCGTATGCCGTCTTCTGCTTG

# Count RAW Reads
```zgrep -c "@M" *.fastq.gz```

* T4-1_S1_L001_R1_001.fastq.gz:3438991
* T4-1_S1_L001_R2_001.fastq.gz:3438991
* T4-6_S1_L001_R1_001.fastq.gz:4070019
* T4-6_S1_L001_R2_001.fastq.gz:4070019
* T4-8_S1_L001_R1_001.fastq.gz:5059437
* T4-8_S1_L001_R2_001.fastq.gz:5059437
* T5-1_S2_L001_R1_001.fastq.gz:3278023
* T5-1_S2_L001_R2_001.fastq.gz:3278023
* T5-6_S2_L001_R1_001.fastq.gz:4508425
* T5-6_S2_L001_R2_001.fastq.gz:4508425
* T5-8_S2_L001_R1_001.fastq.gz:5257196
* T5-8_S2_L001_R2_001.fastq.gz:5257196
* T7-1_S3_L001_R1_001.fastq.gz:3480865
* T7-1_S3_L001_R2_001.fastq.gz:3480865
* T7-6_S3_L001_R1_001.fastq.gz:4179337
* T7-6_S3_L001_R2_001.fastq.gz:4179337
* T7-8_S3_L001_R1_001.fastq.gz:4055930
* T7-8_S3_L001_R2_001.fastq.gz:4055930

#Run FASTQC to examine data quality
```/home/hputnam/programs/FastQC/fastqc /home/hputnam/Mcap_Spawn/Data/Raw_Data/*fastq.gz -o /home/hputnam/Mcap_Spawn/Data/Raw_QC_Files```

#Examine FASTQC Results of raw files
```scp -r hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Data/Raw_QC_Files //Users/hputnam/MyProjects/Montipora_Spawn_Timing/SEQ_Data/Raw_QC_Files```


#Run multicq from raw qc folder to combine results from all files
http://multiqc.info/

```multiqc .```


#Trim Adapters and poor quality
```mkdir /home/hputnam/Mcap_Spawn/Data/cleaned/```

### Used FastqMcf fastq-mcf sequence quality filter, clipping and processor to trim adapters.
https://expressionanalysis.github.io/ea-utils/
https://github.com/ExpressionAnalysis/ea-utils
-o = output
-l = minimum remaining sequence length (=150)
-q = quality threshold causing base removal (=20)
-w = window size for quality trimming (=5)
-x = 'N' bad read percentage causing cycle removal (=10)
-u = Force disable/enable Illumina PF filtering default = auto
-P = Phred-scale default = auto

#FastqMcf 
## Time 4 
```/home/hputnam/programs/ea-utils.1.1.2-806/fastq-mcf -l 16 -q 20 -w 5 -x 10 -u -P 33 \
/home/hputnam/Mcap_Spawn/Refs/barcodes.fa \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T4-1_S1_L001_R1_001.fastq.gz \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T4-1_S1_L001_R2_001.fastq.gz \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T4-1_S1_L001_R1_001_cleaned.fastq \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T4-1_S1_L001_R2_001_cleaned.fastq &>T4-1.log```

```/home/hputnam/programs/ea-utils.1.1.2-806/fastq-mcf -l 16 -q 20 -w 5 -x 10 -u -P 33 \
/home/hputnam/Mcap_Spawn/Refs/barcodes.fa \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T4-6_S1_L001_R1_001.fastq.gz \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T4-6_S1_L001_R2_001.fastq.gz \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T4-6_S1_L001_R1_001_cleaned.fastq \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T4-6_S1_L001_R2_001_cleaned.fastq &>T4-6.log```

```/home/hputnam/programs/ea-utils.1.1.2-806/fastq-mcf -l 16 -q 20 -w 5 -x 10 -u -P 33 \
/home/hputnam/Mcap_Spawn/Refs/barcodes.fa \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T4-8_S1_L001_R1_001.fastq.gz \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T4-8_S1_L001_R2_001.fastq.gz \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T4-8_S1_L001_R1_001_cleaned.fastq \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T4-8_S1_L001_R2_001_cleaned.fastq &>T4-8.log```

## Time 5
```/home/hputnam/programs/ea-utils.1.1.2-806/fastq-mcf -l 16 -q 20 -w 5 -x 10 -u -P 33 \
/home/hputnam/Mcap_Spawn/Refs/barcodes.fa \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T5-1_S2_L001_R1_001.fastq.gz \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T5-1_S2_L001_R2_001.fastq.gz \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T5-1_S2_L001_R1_001_cleaned.fastq \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T5-1_S2_L001_R2_001_cleaned.fastq &>T5-1.log```

```/home/hputnam/programs/ea-utils.1.1.2-806/fastq-mcf -l 16 -q 20 -w 5 -x 10 -u -P 33 \
/home/hputnam/Mcap_Spawn/Refs/barcodes.fa \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T5-6_S2_L001_R1_001.fastq.gz \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T5-6_S2_L001_R2_001.fastq.gz \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T5-6_S2_L001_R1_001_cleaned.fastq \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T5-6_S2_L001_R2_001_cleaned.fastq &>T5-6.log```

```/home/hputnam/programs/ea-utils.1.1.2-806/fastq-mcf -l 16 -q 20 -w 5 -x 10 -u -P 33 \
/home/hputnam/Mcap_Spawn/Refs/barcodes.fa \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T5-8_S2_L001_R1_001.fastq.gz \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T5-8_S2_L001_R2_001.fastq.gz \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T5-8_S2_L001_R1_001_cleaned.fastq \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T5-8_S2_L001_R2_001_cleaned.fastq &>T5-8.log```


## Time 7
```/home/hputnam/programs/ea-utils.1.1.2-806/fastq-mcf -l 16 -q 20 -w 5 -x 10 -u -P 33 \
/home/hputnam/Mcap_Spawn/Refs/barcodes.fa \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T7-1_S3_L001_R1_001.fastq.gz \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T7-1_S3_L001_R2_001.fastq.gz \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T7-1_S3_L001_R1_001_cleaned.fastq \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T7-1_S3_L001_R2_001_cleaned.fastq &>T7-1.log```

```/home/hputnam/programs/ea-utils.1.1.2-806/fastq-mcf -l 16 -q 20 -w 5 -x 10 -u -P 33 \
/home/hputnam/Mcap_Spawn/Refs/barcodes.fa \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T7-6_S3_L001_R1_001.fastq.gz \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T7-6_S3_L001_R2_001.fastq.gz \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T7-6_S3_L001_R1_001_cleaned.fastq \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T7-6_S3_L001_R2_001_cleaned.fastq &>T7-6.log```

```/home/hputnam/programs/ea-utils.1.1.2-806/fastq-mcf -l 16 -q 20 -w 5 -x 10 -u -P 33 \
/home/hputnam/Mcap_Spawn/Refs/barcodes.fa \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T7-8_S3_L001_R1_001.fastq.gz \
/home/hputnam/Mcap_Spawn/Data/Raw_Data/T7-8_S3_L001_R2_001.fastq.gz \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T7-8_S3_L001_R1_001_cleaned.fastq \
-o /home/hputnam/Mcap_Spawn/Data/cleaned/T7-8_S3_L001_R2_001_cleaned.fastq &>T7-8.log```


# Count Cleaned Reads
```grep -c "@M" /home/hputnam/Mcap_Spawn/Data/cleaned/*.fastq```

* /home/hputnam/Mcap_Spawn/Data/cleaned/T4-1_S1_L001_R1_001_cleaned.fastq:3434704
* /home/hputnam/Mcap_Spawn/Data/cleaned/T4-1_S1_L001_R2_001_cleaned.fastq:3434704
* /home/hputnam/Mcap_Spawn/Data/cleaned/T4-6_S1_L001_R1_001_cleaned.fastq:4065700
* /home/hputnam/Mcap_Spawn/Data/cleaned/T4-6_S1_L001_R2_001_cleaned.fastq:4065700
* /home/hputnam/Mcap_Spawn/Data/cleaned/T4-8_S1_L001_R1_001_cleaned.fastq:5053922
* /home/hputnam/Mcap_Spawn/Data/cleaned/T4-8_S1_L001_R2_001_cleaned.fastq:5053922
* /home/hputnam/Mcap_Spawn/Data/cleaned/T5-1_S2_L001_R1_001_cleaned.fastq:3274961
* /home/hputnam/Mcap_Spawn/Data/cleaned/T5-1_S2_L001_R2_001_cleaned.fastq:3274961
* /home/hputnam/Mcap_Spawn/Data/cleaned/T5-6_S2_L001_R1_001_cleaned.fastq:4504000
* /home/hputnam/Mcap_Spawn/Data/cleaned/T5-6_S2_L001_R2_001_cleaned.fastq:4504000
* /home/hputnam/Mcap_Spawn/Data/cleaned/T5-8_S2_L001_R1_001_cleaned.fastq:5251537
* /home/hputnam/Mcap_Spawn/Data/cleaned/T5-8_S2_L001_R2_001_cleaned.fastq:5251537
* /home/hputnam/Mcap_Spawn/Data/cleaned/T7-1_S3_L001_R1_001_cleaned.fastq:3475482
* /home/hputnam/Mcap_Spawn/Data/cleaned/T7-1_S3_L001_R2_001_cleaned.fastq:3475482
* /home/hputnam/Mcap_Spawn/Data/cleaned/T7-6_S3_L001_R1_001_cleaned.fastq:4174182
* /home/hputnam/Mcap_Spawn/Data/cleaned/T7-6_S3_L001_R2_001_cleaned.fastq:4174182
* /home/hputnam/Mcap_Spawn/Data/cleaned/T7-8_S3_L001_R1_001_cleaned.fastq:4051635
* /home/hputnam/Mcap_Spawn/Data/cleaned/T7-8_S3_L001_R2_001_cleaned.fastq:4051635

#Run Fastqc on cleaned files
```mkdir /home/hputnam/Mcap_Spawn/Data/cleaned_QC_Files```

```/home/hputnam/programs/FastQC/fastqc /home/hputnam/Mcap_Spawn/Data/cleaned/*.fastq -o /home/hputnam/Mcap_Spawn/Data/cleaned_QC_Files```

#Examine FASTQC Results of cleaned files
```scp -r hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Data/cleaned_QC_Files /Users/hputnam/MyProjects/Montipora_Spawn_Timing/SEQ_Data/Cleaned_QC_Files```

#Concatenate all R1 and all R2
```cat /home/hputnam/Mcap_Spawn/Data/cleaned/*R1_001_cleaned.fastq  > all_R1_clean.fastq```

```cat /home/hputnam/Mcap_Spawn/Data/cleaned/*R2_001_cleaned.fastq  > all_R2_clean.fastq```

# Count Reads
```grep -c "@M" /home/hputnam/Mcap_Spawn/Data/cleaned/*clean.fastq```

* /home/hputnam/Mcap_Spawn/Data/cleaned/all_R1_clean.fastq:37286123
* /home/hputnam/Mcap_Spawn/Data/cleaned/all_R2_clean.fastq:37286123


#Run Trinity de novo assembly
```/home/hputnam/programs/trinityrnaseq-2.2.0/Trinity --seqType fq  --left /home/hputnam/Mcap_Spawn/Data/cleaned/all_R1_clean.fastq --right /home/hputnam/Mcap_Spawn/Data/cleaned/all_R2_clean.fastq --CPU 40 --max_memory 20G  --min_contig_length 200``` 

# Check Trinity Assembly Stats
```/usr/local/opt/trinityrnaseq/util/TrinityStats.pl /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta > /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.Summary.txt```

#Run BSUSCO on Trinity.fasta 
* http://busco.ezlab.org/files/BUSCO_userguide.pdf

## Assessing assembly and annotation completeness with Benchmarking Universal Single-Copy Orthologs

```python /home/hputnam/programs/BUSCO_v1.22/BUSCO_v1.22.py -o Mcap_Spawn_BUSCO -in /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta -l /home/hputnam/programs/BUSCO_v1.22/eukaryota -m trans```

# Summarizing Trinity Assembly Stats
```/usr/local/opt/trinityrnaseq/util/TrinityStats.pl /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta > Trinity.Sumary.txt```

# Assessing the Read Content of the Transcriptome Assembly

```bowtie2-build /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta```

```bowtie2 -p20 --local --no-unal -x /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta -q -1 /home/hputnam/Mcap_Spawn/Data/cleaned/all_R1_clean.fastq -2 /home/hputnam/Mcap_Spawn/Data/cleaned/all_R2_clean.fastq | samtools view -Sb - | samtools sort -no - - > bowtie2.nameSorted.bam```


# Estimating Transcript Abundance with RSEM in Trinity
 
## Time T4
```/usr/local/opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta --seqType fq --prep_reference --left /home/hputnam/Mcap_Spawn/Data/cleaned/T4-1_S1_L001_R1_001_cleaned.fastq --right /home/hputnam/Mcap_Spawn/Data/cleaned/T4-1_S1_L001_R2_001_cleaned.fastq  --est_method RSEM --aln_method bowtie --trinity_mode --output_prefix T4-1```

```/usr/local/opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta --seqType fq --prep_reference --left /home/hputnam/Mcap_Spawn/Data/cleaned/T4-6_S1_L001_R1_001_cleaned.fastq --right /home/hputnam/Mcap_Spawn/Data/cleaned/T4-6_S1_L001_R2_001_cleaned.fastq  --est_method RSEM --aln_method bowtie --trinity_mode --output_prefix T4-6```

```/usr/local/opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta --seqType fq --prep_reference --left /home/hputnam/Mcap_Spawn/Data/cleaned/T4-8_S1_L001_R1_001_cleaned.fastq --right /home/hputnam/Mcap_Spawn/Data/cleaned/T4-8_S1_L001_R2_001_cleaned.fastq  --est_method RSEM --aln_method bowtie --trinity_mode --output_prefix T4-8```


## Time T5
```/usr/local/opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta --seqType fq --prep_reference --left /home/hputnam/Mcap_Spawn/Data/cleaned/T5-1_S2_L001_R1_001_cleaned.fastq --right /home/hputnam/Mcap_Spawn/Data/cleaned/T5-1_S2_L001_R2_001_cleaned.fastq  --est_method RSEM --aln_method bowtie --trinity_mode --output_prefix T5-1```

```/usr/local/opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta --seqType fq --prep_reference --left /home/hputnam/Mcap_Spawn/Data/cleaned/T5-6_S2_L001_R1_001_cleaned.fastq --right /home/hputnam/Mcap_Spawn/Data/cleaned/T5-6_S2_L001_R2_001_cleaned.fastq  --est_method RSEM --aln_method bowtie --trinity_mode --output_prefix T5-6```

```/usr/local/opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta --seqType fq --prep_reference --left /home/hputnam/Mcap_Spawn/Data/cleaned/T5-8_S2_L001_R1_001_cleaned.fastq --right /home/hputnam/Mcap_Spawn/Data/cleaned/T5-8_S2_L001_R2_001_cleaned.fastq  --est_method RSEM --aln_method bowtie --trinity_mode --output_prefix T5-8```

## Time T7
```/usr/local/opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta --seqType fq --prep_reference --left /home/hputnam/Mcap_Spawn/Data/cleaned/T7-1_S3_L001_R1_001_cleaned.fastq --right /home/hputnam/Mcap_Spawn/Data/cleaned/T7-1_S3_L001_R2_001_cleaned.fastq  --est_method RSEM --aln_method bowtie --trinity_mode --output_prefix T7-1```

```/usr/local/opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta --seqType fq --prep_reference --left /home/hputnam/Mcap_Spawn/Data/cleaned/T7-6_S3_L001_R1_001_cleaned.fastq --right /home/hputnam/Mcap_Spawn/Data/cleaned/T7-6_S3_L001_R2_001_cleaned.fastq  --est_method RSEM --aln_method bowtie --trinity_mode --output_prefix T7-6```

```/usr/local/opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta --seqType fq --prep_reference --left /home/hputnam/Mcap_Spawn/Data/cleaned/T7-8_S3_L001_R1_001_cleaned.fastq --right /home/hputnam/Mcap_Spawn/Data/cleaned/T7-8_S3_L001_R2_001_cleaned.fastq  --est_method RSEM --aln_method bowtie --trinity_mode --output_prefix T7-8```


# Build Transcript and Gene Expression Matrices

/usr/local/opt/trinityrnaseq/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix isoforms_counts_matrix T4-1.isoforms.results T4-6.isoforms.results T4-8.isoforms.results T5-1.isoforms.results T5-6.isoforms.results T5-8.isoforms.results T7-1.isoforms.results T7-6.isoforms.results T7-8.isoforms.results

# Run Differential Expression Analysis

/usr/local/opt/trinityrnaseq/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix /home/hputnam/Mcap_Spawn/RSEM/isoforms_counts_matrix.counts.matrix --method edgeR --samples_file /home/hputnam/Mcap_Spawn/Refs/sample_description.txt 

# Cluster DEG
/usr/local/opt/trinityrnaseq/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix /home/hputnam/Mcap_Spawn/RSEM/isoforms_counts_matrix.TMM.fpkm.matrix -P 0.05 -C 0 --samples /home/hputnam/Mcap_Spawn/Refs/sample_description.txt 


