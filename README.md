# hemoMIPs

The hemoMIPs pipeline is a fast and efficient analysis pipeline for the analysis of multiplexed and targeted NGS datasets created from Molecular Inversion Probes (MIPs). It runs highly automated using conda und snakemake and can be set to use GATK v4 or GATK v3 for variant calling. It reports benign and likely pathogenic variants in a userfriendly HTML report that shows detailed performance statistics and results.

## Pre-requirements

### Conda

The pipeline depends on [Snakemake](https://snakemake.readthedocs.io/en/stable/), a workflow management system that wraps up all scripts and runs them highly automated, in various environments (workstations, clusters, grid, or cloud). Further, we use Conda as software/dependency managment tool. Conda can install snakemake and all neccessary software with its dependencies automatically. Conda installation guidlines can be found here:

https://conda.io/projects/conda/en/latest/user-guide/install/index.html

### Snakemake

After installing Conda, you install Snakemake using Conda and the `environment.yaml` provided in this repository. For this purpose, please clone or download and uncompress the repository first. Then change into the root folder of the local repository. 

```bash
git clone https://github.com/kircherlab/hemoMIPs
cd hemoMIPs
```

We will now initiate three Conda environments which we will need for some preparations as well as getting the Snakemake workflow invoked. The first environment (`hemoMIPs`) will contain only snakemake, the second (`ensemblVEP`) contains Ensembl VEP and htslib, the third (`prepTools`) contains some basic tools for preparing annotations (e.g. bedtools, samtools, htslib, bwa, picard):

```bash
conda env create -n hemoMIPs --file environment.yaml
conda env create -n ensemblVEP --file envs/vep.yml
conda env create -n prepTools --file envs/prep.yml
```

The `ensemblVEP` and `prepTools` environments are only needed for the initial setup and can be deleted afterwards.

### Annotations of Ensembl VEP

We use [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html) to predict variant effects. You will need to install the annotation caches for VEP before you can run the snakmake workflow. For this purpose, you will need to run a tool from the `ensemblVEP` environment that we created above. Please adjust the path to your location in the command line below (`-c vep_cache/`). 

Note that `snakemake` will later install a separate instance of VEP for running the pipeline and that we are only using the above environment to install the caches. Also note, that due to version conflicts with other software, VEP is not included in environments with other software.  If you already have the VEP database, simply adjust the path to your database in the config.yml. We run the pipeline using VEP v98. If you wish to use another version or cache, you should up- or downgrade your specific version of VEP and make sure that the other VEP version is correctly referenced in the workflow.
  
The following commands will download the human VEP cache (approx. 14G), which may take a while. \

```bash
conda activate ensemblVEP
mkdir vep_cache
vep_install -n -a cf -s homo_sapiens -y GRCh37 -c vep_cache/ –CONVERT
conda deactivate
```

## Other required genome annotations

In the following steps, we are preparing the alignment and variant calling index of the reference genome as well as a VCF with known variants. We are using the above created `prepTools` environment:

```bash
conda activate prepTools
```

For alignment, we use the 1000 Genomes phase 2 build of the human reference `hs37d5.fa.gz`, which includes decoy sequences for sequences missing from the assembly. We will need the bwa index and picard/GATK dictionary index of this file.

```bash
mkdir reference_index
cd reference_index
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
bwa index hs37d5.fa
samtools faidx hs37d5.fa
picard CreateSequenceDictionary R=hs37d5.fa O=hs37d5.dict

```

HemoMIPs uses known variants reported by the 1000 Genomes project. To extract known variants for your target region, run the following command using your `target_coords.bed`. Here, we are using the file from the example project:

```bash
cd /~PathTo~/hemoMIPs/
mkdir known_variants
cd known_variants
tabix -h ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20110521/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz -R <( awk 'BEGIN{OFS="\t"}{print $1,$2-50,$3+50}' ../input/example_dataset/target_coords.bed | sort -k1,1 -k2,2n -k3,3n | bedtools merge ) | bgzip -c > phase1_release_v3.20101123.snps_indels_svs.on_target.vcf.gz

tabix -p vcf phase1_release_v3.20101123.snps_indels_svs.on_target.vcf.gz
```

If you decided to include MIPs to capture specific inversion alleles you will also need to provide a BWA index for the inversion MIPs as well as the logic of evaluating those in `scripts/processing/summary_report.py` (lines 138-169). If you do not have inversion probes in your design, set the respective parameter (`Inv`) in the config file to "no". In the following, we will assume that you are using the pipeline for the analysis of the hemophilia MIP design and provide the relevant files with your input folder.

The environments needed to prepare the workflow can be removed at this step. Snakemake will install packages required for the workflow  automatically during the first run of the pipeline. Do not remove the `hemoMIPs` environment as this is needed to invoke snakemake.

```bash
conda deactivate
conda env remove --name ensemblVEP
conda env remove --name prepTools
```

## Config

Almost ready to go. After you prepared files as above, you may need to adjust the locations and names of these files in the `config.yml`. Further, you need to specify your run type, i.e. whether you want to analyze paired-end read or single-end read data as well as your index design (single or double index) in the `config.yml`. The original workflow was developed for paired-end 2 x 120bp with one sample index read. The workflow however allows the analysis of single-ended reads and up to two index reads/technical reads. 

Set the parameters in the config file accordingly: 

```
parameters:
   inv: "yes" #set to "no" when no inversion design is provided
   paired_end_reads: "yes" #set to "no" when single-read sequencing is applied
   double_index: "no" #set to yes when double indexing is applied or a second technical read available 
```

Please note that the workflow supports double indexing, i.e. sequence combinations between the two technical reads identify a specific sample, or the provision of a second technical read (e.g. for unique molecular identifiers, UMIs) which is not used for sample assignment but propogated in a separate BAM field for each read. If you are using double indexing set `double_index` to "yes" and provide a three column `sample_index.lst` file (see below). If sequence information from a technical read should be included, also set `double_index` to "yes", but provide a two column `sample_index.lst` file (see below). In this case, the first index read will be used to assign samples, and the sequence of the second read will be included in the BAM files, but will not be evaluated.

An example config can be found in `example_config.yml`. If you would like to run the example data set, please copy it to `config.yml`:

```bash
cd /~PathTo~/hemoMIPs/
cp example_config.yml config.yml
```

## List of required input files

You need your NGS fastq files together with information about your MIP design and the targeted regions. An example dataset is available with all relevant files in `input/example_dataset`. The required fastq input files should be created using the Illumina `bcl2fastq` tool. The pipeline can handle paired-end and single-end reads with up to two technical reads/index reads (i.e. `Undetermined_S0_L00{lane}_R1_001.fastq.gz`, `Undetermined_S0_L00{lane}_I1_001.fastq.gz`, additional for paired end read data: `Undetermined_S0_L00{lane}_R2_001.fastq.gz`, in case of a second index read: `Undetermined_S0_L00{lane}_I2_001.fastq.gz`). For instance, a paired-end single index dataset could be created by `bcl2fastq --create-fastq-for-index-reads --use-bases-mask 'Y*,I*,Y*'`.

Put your NGS fastq files in input/ together with:
- MIP design file as generated by `https://github.com/shendurelab/MIPGEN` named `hemomips_design.txt`
- Named target regions (coordinates) of your MIP experiment named `target_coords.bed`
- A file containing known benign variants (can be left blank) named `benignVars.txt`. 
- A barcode sample assignment file named `sample_index.lst`

Examples and further information about these files is provided below.

### MIP probe design information

Information about the designed MIP probes and their location in the reference genome is needed as a tab-separated text file for the script `TrimMIParms.py`. The default input file has the following columns: index, score, chr, ext_probe_start, ext_probe_stop, ext_probe_copy, ext_probe_sequence, lig_probe_start, lig_probe_stop, lig_probe_copy, lig_probe_sequence, mip_scan_start_position, mip_scan_stop_position, scan_target_sequence, mip_sequence, feature_start_position, feature_stop_position, probe_strand, failure_flags, gene_name, mip_name. This format is obtained from MIP designs generated by MIPGEN (Boyle et al., 2014), a tool for MIP probe design available on GitHub (https://github.com/shendurelab/MIPGEN). Alternatively, files containing at least the following named columns can be used: chr, ext_probe_start, ext_probe_stop, lig_probe_start, lig_probe_stop, probe_strand, and mip_name. It is critical, that the reported coordinates and chromosome names match the reference genome used in alignment. 

We used Y-chromosome specific targets (SRY) to detect the sex of the patient (see chromosome `Y` in `hemomips_design.txt`). Different Y chromosome targets can be designed for sex determination as the workflow simply counts Y-aligned reads. The pipeline also runs without Y-specific MIPs for sex determination, but in this case will output all samples to be female in the final report.

### Named target regions in BED format

Please describe the target regions of your MIP experiments in a BED file. These regions and names will be used in the HTML report. An example of this BED file is provided below (see also `input/example_dataset/target_coords.bed`):

```
X       154250998       154251277       F8/upstream                             
X       154250827       154250998       F8/5-UTR                                
X       154250674       154250827       F8/1                            
X       154227743       154227906       F8/2                            
...
X       154088696       154088893       F8/25                           
X       154065871       154066037       F8/26                           
X       154064063       154065871       F8/3-UTR                                
X       154064033       154064063       F8/downstream                           
X       138612623       138612894       F9/upstream                             
X       138612894       138612923       F9/5-UTR                                
X       138612923       138613021       F9/1                            
X       138619158       138619342       F9/2                            
...
X       138642889       138643024       F9/7                            
X       138643672       138644230       F9/8                            
X       138644230       138645617       F9/3-UTR                                
X       138645617       138645647       F9/downstream
```

### Known benign variants 

A `benignVars.txt` can be used to describe known benign variants. If no such variants are available, an empty file with this name needs to be provided. If variants are provided in this file, these will be printed in gray font in the HTML report and separated in the CSV output files. An example of the format is provided below. The full file for the hemophila project is available as `input/example_dataset/benignVars.txt`.

```
X_138633280_A/G
X_154065069_T/G
X_138644836_G/A
X_138645058_GT/-
X_138645060_-/GT
X_138645149_T/C
```

### Barcode to sample assignment

A two or three column tab-separated file is required with the sequencing barcode information. The sample name will be used throughout the processing and reporting. If a two colum tab-separated file is provided, the sample barcode sequence is assumed to be in the first index read of the Illumina sequencing run (I1 FastQ read file). The pipeline can also handle double index designs where sequence combinations in the I1 and I2 files identify a specific sample. An example for the sample assignment files is provided below:

Single Index
```
#Seq  Name
ACTGGTAGG	Plate_001_01B.2
GCTCCAACG	Plate_001_01C.3
GCGTAAGAT	Plate_001_01D.4
TGACCATCA	Plate_001_01E.5
GGATTCTCG	Plate_001_01F.6
```

Double Index
```
#Index1	Index2	Name
GGATTCTCG	ACTGGTAGG	Plate_001_01A.1
CATGCGAGA	GCGTAAGAT	Plate_001_01B.2
TGACCATCA	TGACCATCA	Plate_001_01C.3
CATGCGAGA	GGATTCTCG	Plate_001_01D.4
```

## Run pipeline

Ready to go! If you run the pipeline on a cluster see the `cluster.json` for an estimate of minimum requirements for the individual jobs. Note that this depends on your dataset size so you may have to adjust this.
To start the pipeline:

```bash
conda activate hemoMIPs
# dry run to see if everything works
snakemake  --use-conda --configfile config.yml -n
# run the pipeline
snakemake  --use-conda --configfile config.yml
```

We added an example_results folder to enable users to compare the output of the example_dataset analysis to our results.

## Output files

The pipeline outputs varies files in intermediate steps as well as final analysis tables for the user to look at.
Here, we describe the output folder structure. For further information about the various analysis steps, see _Pipeline description_ below.

In the `output/` folder `dataset/` folders (named after your individual datasets) will be generated containing all output files.
Within this folder all processed files can be found in `mapping`, with genotyping files stored in `mapping/gatk4` or `mapping/gatk3`, respectively. The analysis tables and html report files can be found in `report`.

### Mapping
`/output/dataset/mapping/` 

The reads from the primary input fastq files are converted to BAM format (e.g. `mapping/sample.bam`). In case of multiple lanes, these are split into `mapping/sample_lX.bam` files. In these files, overlapping paired-end reads are merged (overlap consensus) and reads are assigned to samples using read groups. Information from the technical reads (I1/I2) is stored in `XI` and `YI` fields for the sequence and `XJ` and `YJ` fields for quality scores. 

Individual (i.e. demultiplexed) sample.bams can be found in `mapping/by_sample/`.

Aligned and MIP arm trimmed files for each sample can be found in BAM format in `mapping/aligned/`. This folder also contains BAM index files. These are index files are for example required to visualize alignments in IGV.

Per sample information about reads aligning to the inversion MIP design (if provided) are stored in `mapping/inversion_mips` as individual BAM files and counts are summarized in `mapping/inversion_mips/inversion_summary_counts.txt`.

#### Genotyping using GATK4
`/output/dataset/mapping/gatk4`

Output files generated by GATK4 HaplotypeCaller (emitting all sites) can be found as `realign_all_samples.bam` and `bam.vcf.gz`. \
Genomic Variant Call Format (GVCF) files for each sample are available in `gatk4/gvcf/` as `SAMPLE.g.vcf.gz` files. \
`realign_all_samples.all_sites.vcf.gz` is the combined VCF generated by GATK4 CombineGVCFs. \
The final genotyped VCF is called `realign_all_samples.vcf.gz`. \
MIP performance statistics can be found in `realign_all_samples.MIPstats.tsv`. \
Variant Effect Pridictions are stored in `realign_all_samples.vep.tsv.gz`.

#### Genotyping using GATK3
`/output/dataset/mapping/gatk3`

This folder contains: \
A realigned BAM generated by GATK3 IndelRealigner: `realign_all_samples.bam`. \
A VCF containing genotypes for all sites generated by GATK3 UnifiedGenotyper: `realign_all_samples.all_sites.vcf.gz`. \
The final VCF with non-homozygote reference alleles: `realign_all_samples.vcf.gz`. \
A filtered list of InDels: `realign_all_samples.indel_check.txt`. \
MIP performance statistics: `realign_all_samples.MIPstats.tsv`. \
Variant Effect Pridictions: `realign_all_samples.vep.tsv.gz`.

## Report
`/output/dataset/report`

Final analysis tables and html files are stored in the `/report/gatk4` or `/report/gatk3` folder depending on which GATK version is used. A description of output files and is available in the sections _Report generation_ and _Report tables in text format_ below.

## Pipeline description

### Primary sequence processing

The primary inputs are raw FastQ files from the sequencing run as well as a sample-to-barcode assignment. In primary processing, reads are converted to BAM format, demultiplexed (storing sample information as read group information), and overlapping paired-end reads are merged and consensus called (Kircher, 2012).

### Alignment and MIP arm trimming

Processed reads are aligned to the reference genome (here GRCh37 build from the 1000 Genomes Project Phase II release) using Burrows-Wheeler Alignment (BWA) 0.7.5 mem (Li and Durbin, 2010). As MIP arm sequence can result in incorrect variant identification (by hiding existing variation below primer sequence), MIP arm sequences are trimmed based on alignment coordinates and new BAM files are created. In this step, we are using MIP design files from MIPgen (Boyle et al., 2014) by default. MIP representation statistics (text output file) are calculated from the aligned files. Further, reads aligning to Y-chromosome-unique probes (SRY) are counted for each sample and reported (text output file). In a separate alignment step, all reads are aligned to a reference sequence file describing only the structural sequence variants as mutant and reference sequences. Results are summarized over all samples with the number of reads aligning to each sequence contig in a text report.

### Coverage analysis and variant calling

Coverage differences between MIPs are handled by down sampling regions of excessive coverage. Variants are genotyped using GATK (McKenna et al., 2010) UnifiedGenotyper (v3.4-46) in combination with IndelRealigner (v3.2-2). Alternatively, GATK v4.0.4.0 HaplotypeCaller is used in gVCF mode in combination with CombineGVCFs and GenotypeGVCFs. The gvcf output files are provided in the `output/dataset/mapping/gatk4/gvcf/` output folder for further sample specific information.

The hemophilia datasets perform similar when run either with the GATK3 or GATK4 workflow. However, in low quality genotype calls the performance might vary and a different call set might be obtained. In a reanalysis performed on one of the hemophilia sequencing experiments, the sample specific genotype agreement is above 0.99 (36 different out of 64,308 genotype calls) between the two GATK versions, with high agreement in associated genotype qualities. We therefore choose GATK4 as the standard setting for the workflow as this versions maintains support, is 50x faster and easier to upgrade.

Variant annotations of the called variants, including variant effect predictions and HGVS variant descriptions are obtained from Ensembl Variant Effect Predictor (McLaren et al., 2016).

### Report generation

Different HTML reports are generated for visualization, interpretation and better access to all information collected in previous steps. There are two entry points to this information, organized as two different HTML reports – one summarizing all variant calls and MIP performance across samples and the other summarizing per-sample results in an overview table. The first report (`summary.html`) provides a more technical sample and variant summary, per region coverage and MIP performance statistics. This report across all samples can be used to assess assay performance (e.g. underperforming MIPs could be redesigned in future assays) and allows identification of suspiciously frequent variants (common variants or systematic errors).

The second report (`report.html`) provides an overview of results for each sample, highlighting putative deleterious variants and taking previously defined common/known benign variants out of focus (gray font). Additional information is provided about potential structural variants and incompletely covered regions. This table also provides an overall sample status field with information about passing and failing samples, as well as flags indicating outlier MIP performances.

Both reports provide links to individual report pages of each sample. The individual reports (`ind_SAMPLENAME.html`), provide quality measures like overall coverage, target region coverage, read counts underlying the inferred sample sex (counting Y aligned reads) and MIP performance statistics (over- or underperforming MIPs in this sample), but most importantly provide detailed information on the identified variants, structural variant call results and regions without coverage (potential deletions). 

### Report tables in text format

In additional to the HTML output files for visualization, results are also presented in computer readable CSV format (comma separated) files. These CSV files can be joined by either the variant or sample specific identifier columns. The following results are summarized in the respective table files:

- `ind_status.csv` outputs the sample sex inferred from SRY counts, reports outlier MIP performance, number of genotype (GT) calls, covered sites within the MIP design regions, average coverage, heterozygous sites, incompletely covered regions, deletions as well as a textual summary in a sample quality flag (e.g. OK, Failed Inversions, Check MIPs). 
- `variant_calls.csv` and `variant_calls_benign.csv` contain all or just benign variants, respectively, with location, genotype, quality scores, allelic depth, coverage and status information. 
- `variant_annotation.csv` provides additional annotations to called variants based on reference and alternative allele information. These annotations include gene name, exonic location, cDNA and CDS position, HGVS Transcript and Protein information, variant rsID, and 1000G allele frequency. 
- `inversion_calls.csv` contains count results for MIPs targeting predefined structural variants. 


## Optional

### GATK v3

GATK v4 is included as a conda environment which automatically installs GATK v4.0.4.0 and all its dependencies.
If you prefer to run the original pipeline using GATK v3 (i.e. GATK 3.2.2 and GATK 3.4-46) you need to change `config.yml` to additionaly include "gatk3" or replace the gatk4 entry. Note that GATK 3.2.2 and 3.4-46 are no longer available for download from the BROAD websites. We therefore provide the required JAR files with this repository rather than obtaining them through Conda. 

### Shed Skin

Shed Skin is an experimental compiler, that can translate pure, but implicitly statically typed Python (2.4-2.6) programs into optimized C++. To fasten (~5x) the read overlapping process one of our python scripts can be translated to C++ with Shed Skin and cross-compiled. This will speed up the analysis but is not crucial for its implementation.

We are providing an example how we were able to cross-compile using shedskin. Please note that this example assumes that miniconda was installed. If you are using another source for Conda, you might need to adjust paths. Further, we need an environment with python v2.6 and the requirements for Shed Skin, which we provide as `envs/shedskin.yml`. Be sure that you are in your root hemoMIPs pipeline folder when executing the following commands.

```bash
# create a new environment
conda env create -f envs/shedskin.yml -n shedskin

mkdir -p ~/miniconda3/envs/shedskin/etc/conda/activate.d
mkdir -p ~/miniconda3/envs/shedskin/etc/conda/deactivate.d

echo '#!/bin/sh
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/shedskin/lib:$LD_LIBRARY_PATH"' > ~/miniconda3/envs/shedskin/etc/conda/activate.d/env_vars.sh

echo '#!/bin/sh
unset LD_LIBRARY_PATH' > ~/miniconda3/envs/shedskin/etc/conda/deactivate.d/env_vars.sh

conda activate shedskin
```
Then download and install Shed Skin v0.9.4 into the bin directory of the hemoMIPs pipeline.

```bash
# Download Shed Skin 0.9.4
wget https://github.com/shedskin/shedskin/releases/download/v0.9.4/shedskin-0.9.4.tgz
# Create bin folder
mkdir -p bin
# Extract Shed Skin and remove file
tar -xzf shedskin-0.9.4.tgz -C bin
rm shedskin-0.9.4.tgz
# Install Shed Skin
cd bin/shedskin-0.9.4
python setup.py install
```
Now we can test the shed Skin installation. We have to modify the `Makefile` to point to the GC library.
```bash
shedskin -L ~/miniconda3/envs/shedskin/include test
sed -i '3s|$| -L ~/miniconda3/envs/shedskin/lib|' Makefile
make
./test
```
The result should look similar to:
```
*** SHED SKIN Python-to-C++ Compiler 0.9.4 ***
Copyright 2005-2011 Mark Dufour; License GNU GPL version 3 (See LICENSE)

[analyzing types..]
********************************100%
[generating c++ code..]
[elapsed time: 1.29 seconds]

hello, world!
```
If the installation or test fails please have a look a the [Shed Skin Dokumentation](https://shedskin.readthedocs.io/en/latest/).

#### Compiling MergeTrimReads.py script

Now we need to compile the `MergeTrimReads.py` script as an extension module using Shed Skin:

```bash
# Go to the script folder
cd scripts/pipeline2.0
# create Makefile and edit it
shedskin -e -L ~/miniconda3/envs/shedskin/include MergeTrimReads
sed -i '3s|$| -L ~/miniconda3/envs/shedskin/lib|' Makefile
# Compile!
make
cd ../../
```

