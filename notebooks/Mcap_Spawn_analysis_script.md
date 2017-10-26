# Analysis for Montipora_Spawn_Timing
Data uploaded and analyzed on Poire 
* Linux poire 3.11.0-26-generic #45~precise1-Ubuntu SMP Tue Jul 15 04:02:35 UTC 2014 x86_64 x86_64 x86_64 GNU/Linux

* @galaxy.geodata.hawaii.edu

### Program List and Versions
* fastqc - v0.11.3
* fastq-mcf - 1.04.803
* multiqc - v1.2
* hmmer - v3.1
* signalp - v4.1
* Trinotate - v3.0.1
* tmhmm - v2.0
* hmmer - v3.1
* ncbi-blast-2.6.0+
* TransDecoder - v3.0.1
* InterProScan - v5.26-65.0
* Panther - v12.0


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

* Add Bundle data

* 1086 and 1088 Putnam 2011

* Putnam et al 2017 (collected 2015)
* https://www.ncbi.nlm.nih.gov/biosample/5607941
* https://www.ncbi.nlm.nih.gov/biosample/5607940


### Used the following adapter seqs to make a MiSeq barcodes file
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

```zgrep -c "@HWUSI" *.fastq.gz```

1086_GTCCGC_L003_R1_001.fastq.gz:1814870
1086_GTCCGC_L003_R2_001.fastq.gz:1814870
1088_GTGAAA_L003_R1_001.fastq.gz:2005551
1088_GTGAAA_L003_R2_001.fastq.gz:2005551

```zgrep -c "@SRR" *.fastq.gz```
SRR4048722.fastq.gz:11945872
SRR4048723.fastq.gz:17931677

### Run FASTQC to examine data quality
```mkdir /home/hputnam/Mcap_Spawn/Data/Raw_QC_Files```

```/home/hputnam/programs/FastQC/fastqc /home/hputnam/Mcap_Spawn/Data/*fastq.gz -o /home/hputnam/Mcap_Spawn/Data/Raw_QC_Files```

```/home/hputnam/programs/FastQC/fastqc /home/hputnam/Mcap_Spawn/Data/Bundles/Raw/*fastq.gz -o /home/hputnam/Mcap_Spawn/Data/Bundles/Raw_QC_Files```

### Examine FASTQC Results of raw files
```scp -r hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Data/Raw_QC_Files /Users/hputnam/MyProjects/Montipora_Spawn_Timing/Data/Raw_QC_File```

```scp -r hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Data/Bundles/Raw_QC_Files /Users/hputnam/MyProjects/Montipora_Spawn_Timing/Data/Raw_QC_File```


### Run multicq from raw qc folder to combine results from all files
http://multiqc.info/

```~/MultiQC/scripts/multiqc .```

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
```mkdir /home/hputnam/Mcap_Spawn/Data/Cleaned_Data```

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


```nohup sh -c 'for file in "1086" "1088"
do
/home/hputnam/programs/ea-utils.1.1.2-806/fastq-mcf \
/home/hputnam/Mcap_Spawn/Refs/barcodes.fa \
/home/hputnam/Mcap_Spawn/Data/Bundles/Raw/${file}_R1.fastq.gz \
/home/hputnam/Mcap_Spawn/Data/Bundles/Raw/${file}_R2.fastq.gz \
-l 50 \
-q 20 \
-w 5 \
-x 10 \
-u \
-P 33 \
-o /home/hputnam/Mcap_Spawn/Data/Bundles/Cleaned/${file}_R1_50_cleaned.fastq.gz \
-o /home/hputnam/Mcap_Spawn/Data/Bundles/Cleaned/${file}_R2_50_cleaned.fastq.gz &> /home/hputnam/Mcap_Spawn/Data/Bundles/Cleaned/${file}_50.log
done'```

```nohup sh -c 'for file in "SRR4048722" "SRR4048723"
do
/home/hputnam/programs/ea-utils.1.1.2-806/fastq-mcf \
/home/hputnam/Mcap_Spawn/Refs/barcodes.fa \
/home/hputnam/Mcap_Spawn/Data/Bundles/Raw/${file}.fastq.gz \
-l 100 \
-q 20 \
-w 5 \
-x 10 \
-u \
-P 33 \
-o /home/hputnam/Mcap_Spawn/Data/Bundles/Cleaned/${file}_100_cleaned.fastq.gz \
&> /home/hputnam/Mcap_Spawn/Data/Bundles/Cleaned/${file}_100.log
done'```

# Count Cleaned Reads
```zgrep -c "@M" /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/*.fastq.gz```
T4-10-includes-adapter_S1_L001_R1_001_cleaned.fastq.gz:4160546
T4-10-includes-adapter_S1_L001_R2_001_cleaned.fastq.gz:4160546
T4-16-includes-adapter_S1_L001_R1_001_cleaned.fastq.gz:5678719
T4-16-includes-adapter_S1_L001_R2_001_cleaned.fastq.gz:5678719
T4-17-includes-adapter_S1_L001_R1_001_cleaned.fastq.gz:5645303
T4-17-includes-adapter_S1_L001_R2_001_cleaned.fastq.gz:5645303
T4-1-includes-adapter_S1_L001_R1_001_cleaned.fastq.gz:3411248
T4-1-includes-adapter_S1_L001_R2_001_cleaned.fastq.gz:3411248
T4-6-includes-adapter_S1_L001_R1_001_cleaned.fastq.gz:4034708
T4-6-includes-adapter_S1_L001_R2_001_cleaned.fastq.gz:4034708
T4-8-includes-adapter_S1_L001_R1_001_cleaned.fastq.gz:5020101
T4-8-includes-adapter_S1_L001_R2_001_cleaned.fastq.gz:5020101
T5-10-includes-adapter_S2_L001_R1_001_cleaned.fastq.gz:5146286
T5-10-includes-adapter_S2_L001_R2_001_cleaned.fastq.gz:5146286
T5-16-includes-adapter_S2_L001_R1_001_cleaned.fastq.gz:5016739
T5-16-includes-adapter_S2_L001_R2_001_cleaned.fastq.gz:5016739
T5-17-includes-adapter_S2_L001_R1_001_cleaned.fastq.gz:4736384
T5-17-includes-adapter_S2_L001_R2_001_cleaned.fastq.gz:4736384
T5-1-includes-adapter_S2_L001_R1_001_cleaned.fastq.gz:3255222
T5-1-includes-adapter_S2_L001_R2_001_cleaned.fastq.gz:3255222
T5-6-includes-adapter_S2_L001_R1_001_cleaned.fastq.gz:4471057
T5-6-includes-adapter_S2_L001_R2_001_cleaned.fastq.gz:4471057
T5-8-includes-adapter_S2_L001_R1_001_cleaned.fastq.gz:5214323
T5-8-includes-adapter_S2_L001_R2_001_cleaned.fastq.gz:5214323
T7-10-includes-adapter_S3_L001_R1_001_cleaned.fastq.gz:3740434
T7-10-includes-adapter_S3_L001_R2_001_cleaned.fastq.gz:3740434
T7-16-includes-adapter_S3_L001_R1_001_cleaned.fastq.gz:4629895
T7-16-includes-adapter_S3_L001_R2_001_cleaned.fastq.gz:4629895
T7-17-includes-adapter_S3_L001_R1_001_cleaned.fastq.gz:4315370
T7-17-includes-adapter_S3_L001_R2_001_cleaned.fastq.gz:4315370
T7-1-includes-adapter_S3_L001_R1_001_cleaned.fastq.gz:3453703
T7-1-includes-adapter_S3_L001_R2_001_cleaned.fastq.gz:3453703
T7-6-includes-adapter_S3_L001_R1_001_cleaned.fastq.gz:4151588
T7-6-includes-adapter_S3_L001_R2_001_cleaned.fastq.gz:4151588
T7-8-includes-adapter_S3_L001_R1_001_cleaned.fastq.gz:4027982
T7-8-includes-adapter_S3_L001_R2_001_cleaned.fastq.gz:4027982

* Reads for assembly
All reads paired end = 160,291,216
Half = 80,109,608

```zgrep -c "@HWUSI" *.fastq.gz```
1086_R1_50_cleaned.fastq.gz:1388212
1086_R2_50_cleaned.fastq.gz:1388212
1088_R1_50_cleaned.fastq.gz:1440311
1088_R2_50_cleaned.fastq.gz:1440311

```zgrep -c "@SRR" SRR*.fastq.gz```
SRR4048722_100_cleaned.fastq.gz:11867209
SRR4048723_100_cleaned.fastq.gz:17828337

#Run Fastqc on cleaned files
```mkdir /home/hputnam/Mcap_Spawn/Data/Cleaned_QC_Files```

```/home/hputnam/programs/FastQC/fastqc /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/*.fastq.gz -o /home/hputnam/Mcap_Spawn/Data/Cleaned_QC_Files```

```/home/hputnam/programs/FastQC/fastqc /home/hputnam/Mcap_Spawn/Data/Bundles/Cleaned/*.fastq.gz -o /home/hputnam/Mcap_Spawn/Data/Bundles/Cleaned_QC_Files```

#Examine FASTQC Results of cleaned files
```scp -r hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Data/Cleaned_QC_Files /Users/hputnam/MyProjects/Montipora_Spawn_Timing/Data```

```scp -r hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Data/Bundles/Cleaned_QC_Files/ /Users/hputnam/MyProjects/Montipora_Spawn_Timing/Data/Cleaned_QC_Files/Cleaned_QC_Files```

```~/MultiQC/scripts/multiqc .```

#Concatenate all R1 and all R2
```cat /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/*R1_001_cleaned.fastq.gz  > all_R1_clean.fastq.gz```

```cat /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/*R2_001_cleaned.fastq.gz  > all_R2_clean.fastq.gz```

# Count Reads
```zgrep -c "@M" /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/*clean.fastq```

read count = 80,109,608

#Run Trinity de novo assembly
* with in silico normalization
```mkdir /home/hputnam/Mcap_Spawn/Assembly```

```nohup /home/hputnam/programs/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq  --left /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/all_R1_clean.fastq.gz --right /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/all_R2_clean.fastq.gz --CPU 40 --max_memory 400G  --min_contig_length 200``` 

* without in silico normalization
```mkdir /home/hputnam/Mcap_Spawn/Assembly/Assembly_no_norm```

```nohup /home/hputnam/programs/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq  --left /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/all_R1_clean.fastq.gz --right /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/all_R2_clean.fastq.gz --CPU 30 --max_memory 20G  --min_contig_length 200 --no_normalize_reads``` 

# Check Trinity Assembly Stats
```/home/hputnam/programs/trinityrnaseq-Trinity-v2.4.0/util/TrinityStats.pl /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta > /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.Summary.txt```

#Run BSUSCO on Trinity.fasta 
* http://busco.ezlab.org/files/BUSCO_userguide.pdf

## Assessing assembly and annotation completeness with Benchmarking Universal Single-Copy Orthologs
```mkdir BUSCO/alv```

```python /home/hputnam/programs/BUSCO_v1.22/BUSCO_v1.22.py -o Mcap_Spawn_BUSCO -in /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta -l /home/hputnam/Mcap_Spawn/Refs/alveolata_stramenophiles_ensembl -m trans```

```mkdir BUSCO/euk```

```python /home/hputnam/programs/BUSCO_v1.22/BUSCO_v1.22.py -o Mcap_Spawn_BUSCO -in /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta -l /home/hputnam/Mcap_Spawn/Refs/eukaryota_odb9 -m trans```

```mkdir BUSCO/meta```

```python /home/hputnam/programs/BUSCO_v1.22/BUSCO_v1.22.py -o Mcap_Spawn_BUSCO -in /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta -l /home/hputnam/Mcap_Spawn/Refs/metazoa_odb9 -m trans```

# Estimating Transcript Abundance 
```mkdir /home/hputnam/Mcap_Spawn/RSEM```

```nohup sh -c 'for file in "T4-10-includes-adapter_S1" "T4-16-includes-adapter_S1" "T4-17-includes-adapter_S1" "T4-1-includes-adapter_S1" "T4-6-includes-adapter_S1" "T4-8-includes-adapter_S1" "T5-10-includes-adapter_S2" "T5-16-includes-adapter_S2" "T5-17-includes-adapter_S2" "T5-1-includes-adapter_S2" "T5-6-includes-adapter_S2" "T5-8-includes-adapter_S2"  "T7-10-includes-adapter_S3" "T7-16-includes-adapter_S3" "T7-17-includes-adapter_S3" "T7-1-includes-adapter_S3" "T7-6-includes-adapter_S3" "T7-8-includes-adapter_S3" 
do
/home/hputnam/programs/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl \
--transcripts /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta \
--seqType fq \
--prep_reference \
--left /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/${file}_L001_R1_001_cleaned.fastq.gz \
--right /home/hputnam/Mcap_Spawn/Data/Cleaned_Data/${file}_L001_R2_001_cleaned.fastq.gz  \
--est_method RSEM \
--aln_method bowtie2 \
--trinity_mode \
--output_dir /home/hputnam/Mcap_Spawn/RSEM \
--output_prefix ${file}
done'```

#### 2011 egg sperm bundles

```nohup sh -c 'for file in "1086" "1088"  
do
/home/hputnam/programs/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl \
--transcripts /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta \
--seqType fq \
--prep_reference \
--left /home/hputnam/Mcap_Spawn/Data/Bundles/Cleaned/${file}_R1_50_cleaned.fastq.gz \
--right /home/hputnam/Mcap_Spawn/Data/Bundles/Cleaned/${file}_R2_50_cleaned.fastq.gz  \
--est_method RSEM \
--aln_method bowtie2 \
--trinity_mode \
--output_dir /home/hputnam/Mcap_Spawn/Data/Bundles/RSEM \
--output_prefix ${file}
done'```

#### 2015 egg sperm bundles

```nohup sh -c 'for file in "SRR4048722" "SRR4048723"  
do
/home/hputnam/programs/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl \
--transcripts /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta \
--seqType fq \
--prep_reference \
--single /home/hputnam/Mcap_Spawn/Data/Bundles/Cleaned/${file}_100_cleaned.fastq.gz \
--est_method RSEM \
--aln_method bowtie2 \
--trinity_mode \
--output_dir /home/hputnam/Mcap_Spawn/Data/Bundles/RSEM \
--output_prefix ${file}
done'```


# Build Transcript and Gene Expression Matrices

```/home/hputnam/programs/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix isoforms_counts_matrix *.isoforms.results```

```/home/hputnam/programs/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method RSEM --out_prefix adult_bundle_isoforms_counts_matrix  /home/hputnam/Mcap_Spawn/RSEM/*.isoforms.results /home/hputnam/Mcap_Spawn/Data/Bundles/RSEM/*.isoforms.results```



# Statistical analysis and plotting with DESeq2
* See R script for statistical analysis with DESeq2

```scp hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/RSEM/isoforms_counts_matrix.counts.matrix /Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Data```

```scp hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Data/Bundles/RSEM/all.isoforms_counts_matrix.counts.matrix /Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Data```

```scp hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Data/Bundles/RSEM/adult_bundle_isoforms_counts_matrix.counts.matrix /Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Data```


# Taxonomic Annotation
# Split into Host and Symbiodinium and other contigs


Coral host reference = Coral.fa
* A_digitifera.fa obtained from OIST Shinzato et al 
* Mcavernosa from Matz lab
* Orbicella from Prada et al 2016
* Stylophora pistillata from Voolstra et al 2017
* Mcapitata

```cat A_digitifera.fa gsd1_racon2.fasta MPSW01.1.fsa_nt Spis.genome.scaffold.final.fa 20170313.mcap.falcon.errd.fasta > Coral.fa```


Symbiodinium reference = Symbiodinium.fa
Pinzon et al and Aranda et al 
* B1 - Shoguchi E et al. 2013, S. minutum (type B1, strain Mf1.05b; 76284 contigs—45263394bp from the host O. faveolata Bayer et al 2012 
* A3 - genomic sequences from cultured Symbiodinium types S. fitti (type A3; 97 259 contigs—21 653 717 bp) Pinzon et al 2015
* C1 - genomic sequences from cultured Symbiodinium types C1 (82 331 contigs—44 078 667 bp) Pinzon et al 2015
* A - transcriptome data from S. microadriaticum (type A, KB8 strain—72152 contigs—61869232bp from the host Cassiopeia spp.) Pinzon et al 2015
* F - Lin et al 2015
* A - Aranda et al 2016

Bacterial reference = Bacteria.fa

```wget "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"```

```awk -F '\t' '{if($12=="Complete Genome") print $20}' assembly_summary.txt > assembly_summary_complete_genomes.txt```


```mkdir GbBac```

```nohup sh -c 'for next in $(cat assembly_summary_complete_genomes.txt); do wget -P GbBac "$next"/*[^m]_genomic.fna.gz; done'```

```gunzip GbBac/*.gz```

```cat GbBac/*.fna > Bacteria.fa```

Viral reference = viruses.fa
http://darwin.informatics.indiana.edu/col/courses/L519/Lab/Lab1/gb2fasta2.pl.txt

```wget "http://www.phantome.org/Downloads/genomes/genbank/2017-06-01.tgz"```

```tar xvzf 2017-06-01.tgz```

```for next in $(ls 2017-06-01); do perl gb2fasta2.pl -gbk 2017-06-01/"$next" -fasta 2017-06-01/"$next".fasta; done```

```cat 2017-06-01/*.fasta > viruses.fa```


##### Make Blast dbs 

```makeblastdb -in Symbiodinium.fa -dbtype nucl```

```makeblastdb -in Coral.fa -dbtype nucl```

```makeblastdb -in Bacteria.fa -dbtype nucl```

```makeblastdb -in viruses.fa -dbtype nucl```

##### BLAST assembly against references

```mkdir Taxa```
```cd Taxa```

```nohup ~/programs/ncbi-blast-2.6.0+/bin/blastn -query ~/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta -db ~/Refs/Symbiodinium.fa -num_threads 30 -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 -outfmt 6 > sym.outfmt6```


```nohup ~/programs/ncbi-blast-2.6.0+/bin/blastn -query ~/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta -db ~/Refs/Coral.fa -num_threads 30 -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 -outfmt 6 > coral.outfmt6```


```nohup ~/programs/ncbi-blast-2.6.0+/bin/blastn -query ~/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta -db ~/Refs/Bacteria.fa -num_threads 30 -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 -outfmt 6 > bacteria.outfmt6```


```nohup ~/programs/ncbi-blast-2.6.0+/bin/blastn -query ~/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta -db ~/Refs/viruses.fa -num_threads 30 -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 -outfmt 6 > viruses.outfmt6```


* header info for blast results
query_id        subject_id      pct_identity    aln_length      n_of_mismatches gap_openings    q_start q_end   s_start   s_end   e_value bit_score

Trinity seqs = 776496

* virus = 54
* bacteria = 1830
* symbiodinium = 38,533
* coral = 555,164
* nohit = 180,915


##### Filter Hits for single assignment per transcript

look for hits that match to multiple databases and create a filter to assign them to only one taxon


awk -F'\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' coral.outfmt6 sym.outfmt6 > coral.sym.matches

awk -F'\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' coral.outfmt6 bacteria.outfmt6 > coral.bac.matches




Explanation

    -F'\t' : sets the field separator to tab.
    NR==FNR : NR is the current input line number and FNR the current file's line number. The two will be equal only while the 1st file is being read.

    c[$1$2]++; next : if this is the 1st file, save the 1st two fields in the c array. Then, skip to the next line so that this is only applied on the 1st file.

    c[$1$2]>0 : the else block will only be executed if this is the second file so we check whether fields 1 and 2 of this file have already been seen (c[$1$2]>0) and if they have been, we print the line. In awk, the default action is to print the line so if c[$1$2]>0 is true, the line will be printed.





# Annotation with Trinotate

#### Build Databases

```cd /home/hputnam/Mcap_Spawn/Annot/Trinotate/```

```/home/hputnam/programs/Trinotate-3.0.1/admin/Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate```

```/home/hputnam/programs/ncbi-blast-2.6.0+/bin/makeblastdb -in uniprot_sprot.pep -dbtype prot```

#### BLASTx

```mkdir BLAST```
```cd BLAST```

```nohup /home/hputnam/programs/ncbi-blast-2.6.0+/bin/blastx -query /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta -db /home/hputnam/Mcap_Spawn/Annot/Trinotate/uniprot_sprot.pep -num_threads 40 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6```

scp hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Annot/blastx.outfmt6 /Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Data```


#### TransDecoder

```mkdir Transdecoder```
```cd Transdecoder```

```nohup /home/hputnam/programs/TransDecoder-3.0.1/TransDecoder.LongOrfs -t /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta```

* Use file: Trinity.fasta.transdecoder_dir/longest_orfs.pep  for Pfam and/or BlastP searches to enable homology-based coding region identification.

* Then, run TransDecoder.Predict for your final coding region predictions.

#### BLASTp
```cd BLAST```

```nohup /home/hputnam/programs/ncbi-blast-2.6.0+/bin/blastp -query /home/hputnam/Mcap_Spawn/Annot/Transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep -db /home/hputnam/Mcap_Spawn/Annot/Trinotate/uniprot_sprot.pep -num_threads 40 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6```

#### HMMER
```mkdir HMMER```
```cd HMMER```

```nohup /home/hputnam/programs/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan --cpu 50 --domtblout TrinotatePFAM.out /home/hputnam/Mcap_Spawn/Annot/Trinotate/Pfam-A.hmm /home/hputnam/Mcap_Spawn/Annot/Transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep > pfam.log```

#### Transdecoder Predict

nohup /home/hputnam/programs/TransDecoder-3.0.1/TransDecoder.Predict --cpu 40 -t /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta  --retain_pfam_hits /home/hputnam/Mcap_Spawn/Annot/HMM TrinotatePFAM.out retain_blastp_hits /home/hputnam/Mcap_Spawn/Annot/BLAST/blastp.outfmt6```

***"The final coding region predictions will now include both those regions that have sequence characteristics consistent with coding regions in addition to those that have demonstrated blast homology or pfam domain content"***


# Use predicted proteins for rest of annotation

#### signalP
```mkdir SignalP```
```cd SignalP```

```nohup /home/hputnam/programs/signalp-4.1/signalp -f short -n signalp.out /home/hputnam/Mcap_Spawn/Annot/Transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep```

* Changed script to max input of 1,000,000 seqs

#### tmHMM
```mkdir tmHMM```
```cd tmHMM```

```nohup /home/hputnam/programs/tmhmm-2.0c/bin/tmhmm --short < /home/hputnam/Mcap_Spawn/Annot/Transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep > tmhmm.out```

#### RNAMMER
```mkdir RNAMMER```
```cd RNAMMER```

```nohup /home/hputnam/programs/Trinotate-3.0.1/util/rnammer_support/RnammerTranscriptome.pl --transcriptome /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta --path_to_rnammer /home/hputnam/programs/rnammer```


# Build Trinotate SQLite Database

* Need to rebuild... 

/home/hputnam/programs/trinityrnaseq-Trinity-v2.4.0/util/support_scripts/get_Trinity_gene_to_trans_map.pl /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta >  Trinity.fasta.gene_trans_map```

#### Load transcripts and coding regions

/home/hputnam/programs/Trinotate-3.0.1/Trinotate Trinotate.sqlite init --gene_trans_map /home/hputnam/Mcap_Spawn/Annot/SQL_DB/Trinity.fasta.gene_trans_map --transcript_fasta /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta --transdecoder_pep /home/hputnam/Mcap_Spawn/Annot/TransDecoder/Trinity.fasta.transdecoder.pep```

#### Load BLAST homologies
/home/hputnam/programs/Trinotate-3.0.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastx /home/hputnam/Mcap_Spawn/Annot/blastx.outfmt6```

/home/hputnam/programs/Trinotate-3.0.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastp /home/hputnam/Mcap_Spawn/Annot/pep/blastp.outfmt6```

#### Load PFAM 

home/hputnam/programs/Trinotate-3.0.1/Trinotate Trinotate.sqlite LOAD_pfam /home/hputnam/Mcap_Spawn/Annot/HMM/TrinotatePFAM.out```


#### Load transmembrane domains

/home/hputnam/programs/Trinotate-3.0.1/Trinotate Trinotate.sqlite LOAD_tmhmm /home/hputnam/Mcap_Spawn/Annot/tmHH/tmhmm.out```

#### Load signal peptide predictions

/home/hputnam/programs/Trinotate-3.0.1/Trinotate Trinotate.sqlite LOAD_signalp /home/hputnam/Mcap_Spawn/Annot/SP/signalp.out```

#### Output Annotation Report

/home/hputnam/programs/Trinotate-3.0.1/Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls```

scp hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Annot/SQL_DB/trinotate_annotation_report.xls /Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Data```

#### Extract GO terms 

/home/hputnam/programs/Trinotate-3.0.1/util/extract_GO_assignments_from_Trinotate_xls.pl \
--Trinotate_xls /home/hputnam/Mcap_Spawn/Annot/SQL_DB/trinotate_annotation_report.xls \
-T \
> /home/hputnam/Mcap_Spawn/Annot/SQL_DB/go_annotations.txt```

sed 's/,/;/g' /home/hputnam/Mcap_Spawn/Annot/SQL_DB/go_annotations.txt > GO_for_MWU.txt```

scp hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Annot/SQL_DB/GO_for_MWU.txt /Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Data```

scp hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Annot/SQL_DB/go_annotations.txt /Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Data```

##### Extract Gene Lengths
/home/hputnam/programs/trinityrnaseq-Trinity-v2.4.0/util/misc/fasta_seq_length.pl  /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta > /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta.seq_lens```

/home/hputnam/programs/trinityrnaseq-Trinity-v2.4.0/util/misc/TPM_weighted_gene_length.py  \
--gene_trans_map /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta.gene_trans_map \
--trans_lengths /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta.seq_lens \
--TPM_matrix /home/hputnam/Mcap_Spawn/RSEM/isoforms_counts_matrix.TMM.EXPR.matrix > /home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.gene_lengths.txt```

scp hputnam@galaxy.geodata.hawaii.edu:/home/hputnam/Mcap_Spawn/Assembly/trinity_out_dir/Trinity.fasta.seq_lens /Users/hputnam/MyProjects/Montipora_Spawn_Timing/RAnalysis/Data```




#### InterProScan
* https://github.com/ebi-pf-team/interproscan
* interproscan-5.23-62.0
* panther-data-12.0 downloaded 21 October 2017
* used data folder from version 5.26-65.0 downloaded 21 October 2017

```tar -xvzf interproscan-5.23-62.0-64-bit.tar.gz```

load panther data

```wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-12.0.tar.gz```

```tar -pxvzf panther-data-12.0.tar.gz```

input file 
* /home/hputnam/Mcap_Spawn/Annot/Transdecoder/Trinity.fasta.transdecoder_dir/longest_orfs.pep

##### IPS Test
```/home/hputnam/programs/my_interproscan/interproscan-5.23-62.0/interproscan.sh -i /home/hputnam/programs/my_interproscan/interproscan-5.23-62.0/test_proteins.fasta -f tsv -dp```

### Prepare data for input to IPS

##### Prepare input file by Removing astrix
```mkdir IPS```

```sed 's_*__g' /home/hputnam/Mcap_Spawn/Annot/TransDecoder/Trinity.fasta.transdecoder.pep > ISP.input.predicted.pep```

##### run a test of IPS on our data 

```nano Mcap_test.pep```

```~/programs/my_interproscan/interproscan-5.23-62.0/interproscan.sh -i Mcap_test.pep -f tsv```

* chose no lookup this -dp disable precalculation call

##### Reset the number of processors used

```nano ~/programs/my_interproscan/interproscan-5.23-62.0/interproscan.properties```

* reset to: 
* Number of embedded workers at start time
number.of.embedded.workers=1
* Maximum number of embedded workers
maxnumber.of.embedded.workers=39

##### Run IPS

```nohup ~/programs/my_interproscan/interproscan-5.23-62.0/interproscan.sh \
-d ~/Mcap_Spawn/Annot/IPS \
-dp \
-i ISP.input.predicted.pep \
-f tsv```


# Integrate Trinotate with IPS











##### Software carpentry SQL query Programming from databases R
* https://swcarpentry.github.io/sql-novice-survey/11-prog-R/
* Load as database
* joins
* tidying
### joins in R



