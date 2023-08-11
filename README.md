Cayman
======

Cayman (<ins>C</ins>arbohydrate <ins>a</ins>ctive enz<ins>y</ins>mes profiling of <ins>m</ins>et<ins>a</ins>ge<ins>n</ins>omes) is a command-line profiling tool for profiling CAZyme abundances in metagenomic datasets. It takes as input (preferably) cleaned -- quality-filtered and host-filtered -- metagenomic shotgun reads and produces a matrix of CAZyme
Reads-Per-Kilobase-Million (RPKM) abundances for your sample. Cayman makes heavy use of the functional profiling library [`gqlib`](https://github.com/cschu/gqlib).

### Prerequisites

#### Dependencies
  - python>=3.7,<3.11
  - bwa=0.7.17
  - numpy==1.24.2
  - pandas==1.5.3
  - pysam
  - gqlib>=2.11.0

  #### Metagenomics reference datasets and CAZyme catalogues

  Cayman uses prevalence-filtered reference data sets from the [Global Microbial Gene Catalog (GMGC)](https://gmgc.embl.de/). We annotated these datasets with our dedicated CAZyme annotation method (cf. Ducarmon & Karcher et al.). The filtered GMGC datasets and their corresponding CAZyme annotations can be downloaded from Zenodo.

  [TABLE]

  Prior to your first profiling run, you will have to build a bwa index from the respective GMGC reference dataset.

  ```
  bwa index -p <index_name> [-b blocksize] /path/to/dataset
  ```

  If you have enough memory available, setting `-b` to a higher value than the default (`10,000,000`), e.g. `100,000,000`, [may speed up the index generation](https://github.com/lh3/bwa/issues/104).



## Installation
Cayman can most easily be installed using ....

<!-- For your biome of interest, you will have to download the respective gene catalog and its CAZyme annotation file, which can be found on Zenodo under the following identifier:  -->


## Running Cayman

Cayman can be run from the command line as follows:

```
cayman \
  <input_options> \
  --annotation_db </path/to/db> \
  --bwa_index </path/to/bwa_index> \
  [--out_prefix <prefix>] \
  [--min_identity <float>] \
  [--min_seqlen <int>] \
  [--cpus_for_alignment <int>]
```

### Mandatory parameters

* `<input_options>`

  1. Read files need to be in fastq format (best with `fastq` or `fq` file ending) and can be gzip compressed.
  2. The `<input_options>` parameters depend on the library layout of your samples:
      * Paired-end data can be specified with `--reads1 </path/to/reads1> --reads2 </path/to/reads2>`.
      * Single-end data can be specified with `--singles </path/to/reads>`.
      * Orphaned reads, i.e. paired-end reads that have lost their mate during an upstream quality control step, can be specified with `--orphans </path/to/orphans>`.

 
  3. Samples comprising multiple fastq files (e.g. from multiple lanes) can be specified as comma-separated lists. In the case of paired-end reads, ensure that the order of the files matches (e.g. `--reads1 sampleX_lane1_R1.fq,sampleX_lane2_R1.fq --reads2 sampleX_lane1_R2.fq,sampleX_lane2_R2.fq`)!


  4. The choice of assigning an unpaired read set to be "true" single-end reads or orphan reads influences the read count distribution.

      * A read pair gets assigned a count of `2 x 0.5 = 1` (as both reads of a pair are derived from the same nucleic acid fragment.)
      * An orphan read gets assigned a count of `1 x 0.5 = 0.5`.
      * An read from a single-end library gets assigned a count of `1`.
  

* `--annotation_db` is the path to a bed4 database containing the reference domain annotation. (format: contig,start,end,domain-type). This contains all the CAZy domain annotations for all ORFs in our gene catalog.

* `--bwa_index` refers to the path of the gene catalog bwa index.

### Optional parameters

* `--out_prefix` here you can define the text string prefix for your output files. If you want to store it in an output folder, then name it /path/to/folder/some_prefix; without some_prefix it will generated hidden files (since they start with a .).

* `--min_identity` refers to the minimum identity level (default 0.97) of the alignment of your read to a CAZyme domain for it to be included.
  
* `--min_seqlen` refers to the minimum length of the alignment to be included (default 45bp).

* `--cpus_for_alignment` the number of cpus to use for alignment.


## Results
- `<out_prefix>_cazy.combined_rpkm.txt` Which contains the sum of both unique and ambigious alignments.
- `<out_prefix>_cazy.unique_rpkm.txt` Which contains the sum of only unique alignments.
- `<out_prefix>_aln_stats.txt` Which contains statistics on the alignment rates

## Reference
- Cayman paper
- A link to gff_quant?
- https://gmgc.embl.de/
