## Overview
Cayman (Carbohydrate active enzymes profiling of metagenomes) is a command-line profiling tool that takes as input cleaned (quality-filtered and host-filtered) metagenomic shotgun reads and gives as output a matrix of CAZyme
Reads-Per-Kilobase-Million (RPKM) abundances for your sample.

## Prerequisites & dependencies
- python>=3.7,<3.11
- bedtools=2.30.0
- bwa=0.7.17
- samtools=1.13
- numpy==1.24.2
- pandas==1.5.3
- pysam
- gqlib>=2.11.0

## Installation
Cayman can most easily be installed using ....

For your biome of interest, you will have to download the respective gene catalog and its CAZyme annotation file, which can be found on Zenodo under the following identifier

## Running Cayman
After installing Cayman, you can run it from the command-line by providing it with sthogun metagenomic reads as follows:

`cayman --reads1 --reads2 --orphans --out_prefix --min_identity --min_seqlen`

Where 
- `--reads1` indicates your forward reads, `--reads2` indicates your reverse reads and `--orphans` indicates orphan reads which may have lost their read pair mate during quality filtering and / or host filtering.
Alternatively, there is also the `--singles` option in case of single-end sequencing. 

- `--min_identity` refers to the minimum identity level (default 0.97) of the alignment of your read to a CAZyme domain for it to be included and `--min_seqlen` refers to the minimum length of the alignment to be included (default 45bp).

- `--out_prefix` here you can define the text string prefix for your output files.

## Results
- `cazy.combined_rpkm.txt` Which contains the sum of both unique and ambigious alignments.
- `cazy.unique_rpkm.txt` Which contains the sum of only unique alignments.
- `aln_stats.txt` Which contains statistics on the alignment rates

## Reference
- Cayman paper
- A link to gff_quant?
- https://gmgc.embl.de/
