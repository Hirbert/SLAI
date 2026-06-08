# SLAI
source code v1.0

# SLAI Installation and Testing Guide

## 1. Download

### 1.1 Source Code

Create a working directory, e.g., `SLAI_test`:

```bash
cd SLAI_test
git clone https://github.com/Hirbert/SLAI.git
cd SLAI
```

### 1.2 UCSC hg19 Reference Sequence

The code requires the UCSC hg19 reference sequence. If this version is not available on your server, download it from the link below and build the corresponding index files. The downloaded `.fa` file should be placed in the `./test/database` directory.

- Download URL: [https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)

## 2. Installation

### 2.1 Installation via Docker

**Build the Docker image:**

```bash
docker build -t slai_test ./
```
**Test the instalation:**
```bash
docker run --rm -v `pwd`/test_data:/data slai_test python3 SLAI.py -h
```
If the following help message appears, the installation is successful:
```text
usage: SLAI.py [-h] --bamfile BAMFILE --fasta FASTA [--outdir OUTDIR]
               --outfile OUTFILE --snv SNV [--flank FLANK] [--cpu CPU]
               [--final FINAL]
```
### 2.2 Direct Installation

**Requirement:** Python version >= 3.6

**Install the required package:**

```bash
pip install pysam==0.19.0
```
**Test the instalation:**
```bash
python3 SLAI.py -h
```
If the following help message appears, the installation is successful:
```text
usage: SLAI.py [-h] --bamfile BAMFILE --fasta FASTA [--outdir OUTDIR]
               --outfile OUTFILE --snv SNV [--flank FLANK] [--cpu CPU]
               [--final FINAL]
```

## 3. Testing with Demo Data

### 3.1 Test Command for Docker Installation

```bash
docker run --rm -v `pwd`/test_data:/data slai python3 SLAI.py \
    --bamfile /data/sample_data/P265700-WatchMaker.demo.bam  \
    --fasta /data/database/hg19.fa \
    --outdir /data/output \
    --snv /data/sample_data/P265700-WatchMaker.snv.chr1.txt \
    --flank 300 \
    --cpu 8 \
    --outfile P265700-WatchMaker.chr1
```
### 3.2 Test Command for Direct Installation

```bash
python3 SLAI.py \
    --bamfile /data/sample_data/P265700-WatchMaker.demo.bam  \
    --fasta HG19_FASTA_FILE_PATH \
    --outdir /data/output \
    --snv /data/sample_data/P265700-WatchMaker.snv.chr1.txt \
    --flank 300 \
    --cpu 8 \
    --outfile P265700-WatchMaker.chr1
```