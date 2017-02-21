# Analysis for Montipora_Spawn_Timing
Data uploaded and analyzed on Poire 
ssh hputnam@galaxy.geodata.hawaii.edu

### Upload data to server
```scp -r /Volumes/NGS_DATA/Hawaii_Mcap/Spawn_2016 hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Data```

### Check upload completion
* match Original 20161208_checksum.md5 with new upload checksum file
 
```md5sum *.gz > 20170220_server_checksum.md5```

>20161208_checksum.md5

* MD5 (T4-1-includes-adapter_S1_L001_R1_001.fastq.gz) = 0d034de745024ef1f3a0a808fb814da2
* MD5 (T4-1-includes-adapter_S1_L001_R2_001.fastq.gz) = e6a8619819ae4ebbda76fab1c16d2ee5
* MD5 (T4-10-includes-adapter_S1_L001_R1_001.fastq.gz) = ee119a26764c90fc2f958d104a75fd9f
* MD5 (T4-10-includes-adapter_S1_L001_R2_001.fastq.gz) = 3c572aa147718915633c86e707269aba
* MD5 (T4-16-includes-adapter_S1_L001_R1_001.fastq.gz) = 23369f48d9c8a09cdcbf3360c5d98af9
* MD5 (T4-16-includes-adapter_S1_L001_R2_001.fastq.gz) = e3a48c6186b0869b94ccc70df65310cd
* MD5 (T4-17-includes-adapter_S1_L001_R1_001.fastq.gz) = 42e199deae9b280dc0c4c46c3eace29b
* MD5 (T4-17-includes-adapter_S1_L001_R2_001.fastq.gz) = c539fd52b0e5835007b81949bec37898
* MD5 (T4-6-includes-adapter_S1_L001_R1_001.fastq.gz) = d2524a6fd255b15b3e3fb2dfe04b8912
* MD5 (T4-6-includes-adapter_S1_L001_R2_001.fastq.gz) = 9eddbe1ab57021043a2971c07c3f9285
* MD5 (T4-8-includes-adapter_S1_L001_R1_001.fastq.gz) = ecb64ffbedac95c4e3cbe0d473ef2236
* MD5 (T4-8-includes-adapter_S1_L001_R2_001.fastq.gz) = 68755f4a785e6b5723674370dfc3974d
* MD5 (T5-1-includes-adapter_S2_L001_R1_001.fastq.gz) = 06eed3df093f468027e006cd61b71c1a
* MD5 (T5-1-includes-adapter_S2_L001_R2_001.fastq.gz) = f462023d1207004cd4d0e1145407bbff
* MD5 (T5-10-includes-adapter_S2_L001_R1_001.fastq.gz) = 82cfe360bbba75e543aa89824e0beab4
* MD5 (T5-10-includes-adapter_S2_L001_R2_001.fastq.gz) = 324d781e5a9897b745311d258871ebf3
* MD5 (T5-16-includes-adapter_S2_L001_R1_001.fastq.gz) = c7b6251c673bad1a5c41358afb6432d8
* MD5 (T5-16-includes-adapter_S2_L001_R2_001.fastq.gz) = 1e97f1615681b512901abc1a791617ae
* MD5 (T5-17-includes-adapter_S2_L001_R1_001.fastq.gz) = e28b444af3380762f8181e49b3ec8639
* MD5 (T5-17-includes-adapter_S2_L001_R2_001.fastq.gz) = a862e594e4c8a2f0ef0901d9f833777b
* MD5 (T5-6-includes-adapter_S2_L001_R1_001.fastq.gz) = 2eeb28b45c5a000577755c830e36b7be
* MD5 (T5-6-includes-adapter_S2_L001_R2_001.fastq.gz) = 020049d0cda29af60e90925afaa3b0c7
* MD5 (T5-8-includes-adapter_S2_L001_R1_001.fastq.gz) = 4e8e1b60f5f8aea8769133e0f95e36c1
* MD5 (T5-8-includes-adapter_S2_L001_R2_001.fastq.gz) = 4e62cb401dadfff0b8ebdd5e79a67348
* MD5 (T7-1-includes-adapter_S3_L001_R1_001.fastq.gz) = 0a1c2963ed23888ee7fa18ce00f52497
* MD5 (T7-1-includes-adapter_S3_L001_R2_001.fastq.gz) = 608079eadce5498ccab50c72d4ebef2c
* MD5 (T7-10-includes-adapter_S3_L001_R1_001.fastq.gz) = 072e8d420b9cfb6939ffde8f5c264f5e
* MD5 (T7-10-includes-adapter_S3_L001_R2_001.fastq.gz) = 7a2e6ae58b0c6e77ced5cb2fae082d4f
* MD5 (T7-16-includes-adapter_S3_L001_R1_001.fastq.gz) = ef4537d0d48b7dd8b27f0abe4cbd17fd
* MD5 (T7-16-includes-adapter_S3_L001_R2_001.fastq.gz) = a3ab869746c91eaacb6d402d47c7523b
* MD5 (T7-17-includes-adapter_S3_L001_R1_001.fastq.gz) = f94c1468d7de02729c98b93f5f38347a
* MD5 (T7-17-includes-adapter_S3_L001_R2_001.fastq.gz) = 3446a006926122041431f708cd6291d8
* MD5 (T7-6-includes-adapter_S3_L001_R1_001.fastq.gz) = 0722ee9f7a97367403667f878329b9ba
* MD5 (T7-6-includes-adapter_S3_L001_R2_001.fastq.gz) = 0202e6a634547b3275d318ebdf20dfe1
* MD5 (T7-8-includes-adapter_S3_L001_R1_001.fastq.gz) = 820397dbcc084832ad34ba177cea411a
* MD5 (T7-8-includes-adapter_S3_L001_R2_001.fastq.gz) = 5c8792e297e619d97a8f65cc1d2e6b24


> 20170220_server_checksum.md5

* ee119a26764c90fc2f958d104a75fd9f  T4-10-includes-adapter_S1_L001_R1_001.fastq.gz
* 3c572aa147718915633c86e707269aba  T4-10-includes-adapter_S1_L001_R2_001.fastq.gz
* 23369f48d9c8a09cdcbf3360c5d98af9  T4-16-includes-adapter_S1_L001_R1_001.fastq.gz
* e3a48c6186b0869b94ccc70df65310cd  T4-16-includes-adapter_S1_L001_R2_001.fastq.gz
* 42e199deae9b280dc0c4c46c3eace29b  T4-17-includes-adapter_S1_L001_R1_001.fastq.gz
* c539fd52b0e5835007b81949bec37898  T4-17-includes-adapter_S1_L001_R2_001.fastq.gz
* 0d034de745024ef1f3a0a808fb814da2  T4-1-includes-adapter_S1_L001_R1_001.fastq.gz
* e6a8619819ae4ebbda76fab1c16d2ee5  T4-1-includes-adapter_S1_L001_R2_001.fastq.gz
* d2524a6fd255b15b3e3fb2dfe04b8912  T4-6-includes-adapter_S1_L001_R1_001.fastq.gz
* 9eddbe1ab57021043a2971c07c3f9285  T4-6-includes-adapter_S1_L001_R2_001.fastq.gz
* ecb64ffbedac95c4e3cbe0d473ef2236  T4-8-includes-adapter_S1_L001_R1_001.fastq.gz
* 68755f4a785e6b5723674370dfc3974d  T4-8-includes-adapter_S1_L001_R2_001.fastq.gz
* 82cfe360bbba75e543aa89824e0beab4  T5-10-includes-adapter_S2_L001_R1_001.fastq.gz
* 324d781e5a9897b745311d258871ebf3  T5-10-includes-adapter_S2_L001_R2_001.fastq.gz
* c7b6251c673bad1a5c41358afb6432d8  T5-16-includes-adapter_S2_L001_R1_001.fastq.gz
* 1e97f1615681b512901abc1a791617ae  T5-16-includes-adapter_S2_L001_R2_001.fastq.gz
* e28b444af3380762f8181e49b3ec8639  T5-17-includes-adapter_S2_L001_R1_001.fastq.gz
* a862e594e4c8a2f0ef0901d9f833777b  T5-17-includes-adapter_S2_L001_R2_001.fastq.gz
* 06eed3df093f468027e006cd61b71c1a  T5-1-includes-adapter_S2_L001_R1_001.fastq.gz
* f462023d1207004cd4d0e1145407bbff  T5-1-includes-adapter_S2_L001_R2_001.fastq.gz
* 2eeb28b45c5a000577755c830e36b7be  T5-6-includes-adapter_S2_L001_R1_001.fastq.gz
* 020049d0cda29af60e90925afaa3b0c7  T5-6-includes-adapter_S2_L001_R2_001.fastq.gz
* 4e8e1b60f5f8aea8769133e0f95e36c1  T5-8-includes-adapter_S2_L001_R1_001.fastq.gz
* 4e62cb401dadfff0b8ebdd5e79a67348  T5-8-includes-adapter_S2_L001_R2_001.fastq.gz
* 072e8d420b9cfb6939ffde8f5c264f5e  T7-10-includes-adapter_S3_L001_R1_001.fastq.gz
* 7a2e6ae58b0c6e77ced5cb2fae082d4f  T7-10-includes-adapter_S3_L001_R2_001.fastq.gz
* ef4537d0d48b7dd8b27f0abe4cbd17fd  T7-16-includes-adapter_S3_L001_R1_001.fastq.gz
* a3ab869746c91eaacb6d402d47c7523b  T7-16-includes-adapter_S3_L001_R2_001.fastq.gz
* f94c1468d7de02729c98b93f5f38347a  T7-17-includes-adapter_S3_L001_R1_001.fastq.gz
* 3446a006926122041431f708cd6291d8  T7-17-includes-adapter_S3_L001_R2_001.fastq.gz
* 0a1c2963ed23888ee7fa18ce00f52497  T7-1-includes-adapter_S3_L001_R1_001.fastq.gz
* 608079eadce5498ccab50c72d4ebef2c  T7-1-includes-adapter_S3_L001_R2_001.fastq.gz
* 0722ee9f7a97367403667f878329b9ba  T7-6-includes-adapter_S3_L001_R1_001.fastq.gz
* 0202e6a634547b3275d318ebdf20dfe1  T7-6-includes-adapter_S3_L001_R2_001.fastq.gz
* 820397dbcc084832ad34ba177cea411a  T7-8-includes-adapter_S3_L001_R1_001.fastq.gz
* 5c8792e297e619d97a8f65cc1d2e6b24  T7-8-includes-adapter_S3_L001_R2_001.fastq.gz

### Used the following adapter files to make barcodes file
>TruSeq_universal_F
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
> Genomic_DNA_oligonucleotide_sequences_Adapters
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
> Genomic_DNA_oligonucleotide_sequences_Adapters_Read2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
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

### Count RAW Reads
```zgrep -c "@M" *.fastq.gz```

T4-10-includes-adapter_S1_L001_R1_001.fastq.gz:4192400
T4-10-includes-adapter_S1_L001_R2_001.fastq.gz:4192400
T4-16-includes-adapter_S1_L001_R1_001.fastq.gz:5728309
T4-16-includes-adapter_S1_L001_R2_001.fastq.gz:5728309
T4-17-includes-adapter_S1_L001_R1_001.fastq.gz:5685333
T4-17-includes-adapter_S1_L001_R2_001.fastq.gz:5685333
T4-1-includes-adapter_S1_L001_R1_001.fastq.gz:3438991
T4-1-includes-adapter_S1_L001_R2_001.fastq.gz:3438991
T4-6-includes-adapter_S1_L001_R1_001.fastq.gz:4070019
T4-6-includes-adapter_S1_L001_R2_001.fastq.gz:4070019
T4-8-includes-adapter_S1_L001_R1_001.fastq.gz:5059437
T4-8-includes-adapter_S1_L001_R2_001.fastq.gz:5059437
T5-10-includes-adapter_S2_L001_R1_001.fastq.gz:5182807
T5-10-includes-adapter_S2_L001_R2_001.fastq.gz:5182807
T5-16-includes-adapter_S2_L001_R1_001.fastq.gz:5054983
T5-16-includes-adapter_S2_L001_R2_001.fastq.gz:5054983
T5-17-includes-adapter_S2_L001_R1_001.fastq.gz:4777150
T5-17-includes-adapter_S2_L001_R2_001.fastq.gz:4777150
T5-1-includes-adapter_S2_L001_R1_001.fastq.gz:3278023
T5-1-includes-adapter_S2_L001_R2_001.fastq.gz:3278023
T5-6-includes-adapter_S2_L001_R1_001.fastq.gz:4508425
T5-6-includes-adapter_S2_L001_R2_001.fastq.gz:4508425
T5-8-includes-adapter_S2_L001_R1_001.fastq.gz:5257196
T5-8-includes-adapter_S2_L001_R2_001.fastq.gz:5257196
T7-10-includes-adapter_S3_L001_R1_001.fastq.gz:3767511
T7-10-includes-adapter_S3_L001_R2_001.fastq.gz:3767511
T7-16-includes-adapter_S3_L001_R1_001.fastq.gz:4661903
T7-16-includes-adapter_S3_L001_R2_001.fastq.gz:4661903
T7-17-includes-adapter_S3_L001_R1_001.fastq.gz:4349113
T7-17-includes-adapter_S3_L001_R2_001.fastq.gz:4349113
T7-1-includes-adapter_S3_L001_R1_001.fastq.gz:3480865
T7-1-includes-adapter_S3_L001_R2_001.fastq.gz:3480865
T7-6-includes-adapter_S3_L001_R1_001.fastq.gz:4179337
T7-6-includes-adapter_S3_L001_R2_001.fastq.gz:4179337
T7-8-includes-adapter_S3_L001_R1_001.fastq.gz:4055930
T7-8-includes-adapter_S3_L001_R2_001.fastq.gz:4055930

### Run FASTQC to examine data quality
```mkdir /home/hputnam/Mcap_Spawn/Data/Raw_QC_Files```

```/home/hputnam/programs/FastQC/fastqc /home/hputnam/Mcap_Spawn/Data/*fastq.gz -o /home/hputnam/Mcap_Spawn/Data/Raw_QC_Files```

### Examine FASTQC Results of raw files
```scp -r hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Data/Raw_QC_Files /Users/hputnam/MyProjects/Montipora_Spawn_Timing/Data/Raw_QC_File```


### Run multicq from raw qc folder to combine results from all files
http://multiqc.info/

```multiqc .```

### Trim Adapters and poor quality
```mkdir /home/hputnam/Mcap_Spawn/Data/cleaned/```

### Used FastqMcf fastq-mcf sequence quality filter, clipping and processor to trim adapters.
https://expressionanalysis.github.io/ea-utils/
https://github.com/ExpressionAnalysis/ea-utils
-o = output
-l = minimum remaining sequence length (=100)
-q = quality threshold causing base removal (=20)
-w = window size for quality trimming (=5)
-x = 'N' bad read percentage causing cycle removal (=10)
-u = Force disable/enable Illumina PF filtering default = auto
-P = Phred-scale default = auto

### FastqMcf 
mkdir /home/hputnam/Mcap_Spawn/Data/Cleaned_Data

```nohup sh -c 'for file in "T4-10-includes-adapter_S1" "T4-16-includes-adapter_S1" "T4-17-includes-adapter_S1" "T4-1-includes-adapter_S1" "T4-6-includes-adapter_S1" "T4-8-includes-adapter_S1" "T5-10-includes-adapter_S2" "T5-16-includes-adapter_S2" "T5-17-includes-adapter_S2" "T5-1-includes-adapter_S2" "T5-6-includes-adapter_S2" "T5-8-includes-adapter_S2"  "T7-10-includes-adapter_S3" "T7-16-includes-adapter_S3" "T7-17-includes-adapter_S3" "T7-1-includes-adapter_S3" "T7-6-includes-adapter_S3" "T7-8-includes-adapter_S3" 
do
/home/hputnam/programs/ea-utils.1.1.2-806/fastq-mcf \
/home/hputnam/Mcap_Spawn/Refs/barcodes.fa \
/home/hputnam/Mcap_Spawn/Data/${file}_L001_R1_001.fastq.gz \
/home/hputnam/Mcap_Spawn/Data/${file}_L001_R2_001.fastq.gz \
-l 100 \
-q 20 \
-w 5 \
-x 10 \
-u \
-P 33 \
-o /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/${file}_L001_R1_001_cleaned.fastq.gz \
-o /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/${file}_L001_R2_001_cleaned.fastq.gz &> /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/${file}.log
done'```

# Count Cleaned Reads
```grep -c "@M" /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/*.fastq```

#Run Fastqc on cleaned files
```mkdir /home/hputnam/Mcap_Spawn/Data/Cleaned_QC_Files```

```/home/hputnam/programs/FastQC/fastqc /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/*.fastq -o /home/hputnam/Mcap_Spawn/Data/Cleaned_QC_File```

#Examine FASTQC Results of cleaned files
```scp -r hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Data/Cleaned_QC_File /Users/hputnam/MyProjects/Montipora_Spawn_Timing/Data```

#Concatenate all R1 and all R2
```cat /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/*R1_001_cleaned.fastq.gz  > all_R1_clean.fastq```

```cat /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/*R2_001_cleaned.fastq.gz  > all_R1_clean.fastq```

# Count Reads
```grep -c "@M" /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/*clean.fastq```


#Run Trinity de novo assembly
```~/programs/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq  --left /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/all_R1_clean.fastq --right /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/all_R2_clean.fastq --CPU 30 --max_memory 20G  --min_contig_length 200``` 

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


# Estimating Transcript Abundance 
need to add a loop 
## Time T4
```/usr/local/opt/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts /home/hputnam/Mcap_Spawn/trinity_out_dir/Trinity.fasta --seqType fq --prep_reference --left /home/hputnam/Mcap_Spawn/Data/cleaned/T4-1_S1_L001_R1_001_cleaned.fastq --right /home/hputnam/Mcap_Spawn/Data/cleaned/T4-1_S1_L001_R2_001_cleaned.fastq  --est_method RSEM --aln_method bowtie --trinity_mode --output_prefix T4-1```


# Build Transcript and Gene Expression Matrices

/usr/local/opt/trinityrnaseq/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix isoforms_counts_matrix T4-1.isoforms.results T4-6.isoforms.results T4-8.isoforms.results T5-1.isoforms.results T5-6.isoforms.results T5-8.isoforms.results T7-1.isoforms.results T7-6.isoforms.results T7-8.isoforms.results

# Run Differential Expression Analysis

/usr/local/opt/trinityrnaseq/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix /home/hputnam/Mcap_Spawn/RSEM/isoforms_counts_matrix.counts.matrix --method edgeR --samples_file /home/hputnam/Mcap_Spawn/Refs/sample_description.txt 

# Cluster DEG_fpkm
/usr/local/opt/trinityrnaseq/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix /home/hputnam/Mcap_Spawn/RSEM/isoforms_counts_matrix.TMM.fpkm.matrix -P 0.05 -C 0 --samples /home/hputnam/Mcap_Spawn/Refs/sample_description.txt 


# Cluster Expression Profiles
/usr/local/opt/trinityrnaseq/Analysis/DifferentialExpression//define_clusters_by_cutting_tree.pl \
-R  diffExpr.P0.05_C0.matrix.RData --Ptree 60

