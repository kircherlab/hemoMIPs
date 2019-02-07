# hemoMIPs

The hemoMIPs pipeline is a fast and efficient analysis pipeline for MIP targeted NGS-datasets. It runs highly automated using conda und snakemake and can be set to use GATK v3 or GATK v4 for data processing.

## Pre-requirements
=======

module load python/2.7.3
module load pysam/0.7.5
module load gmp/5.0.2 mpfr/3.1.0 mpc/0.8.2 gcc/4.9.1


### Conda

The pipeline depends on snakemake (that wraps up the scripts and runs them highly automated). Snakemake relies on conda to install snakemake itself and most of the depended software needed. Conda installation guidlines can be found here:

https://conda.io/docs/user-guide/install/index.html

### Snakemake

Next, you just have to install Snakemake using Conda. This can be done by creating the main environment to run the hemoMIPs pipeline via Snakemake:

```bash
conda env create -n hemoMIPs --file environment.yaml
```

Finally you can activate the environment via `source activate hemoMIPs`. Now the `snakemake` command is available.

### Shed Skin

Shed Skin is an experimental compiler, that can translate pure, but implicitly statically typed Python (2.4-2.6) programs into optimized C++. To fasten the read overlapping process one of our python scripts have to be translated to C++ with Shed Skin.
This will speed up the analysis drastically but is not crucial for the implementation.
First we need an environment with python v2.6 and the requirements for Shed Skin. Therefore we created the environment file `envs/shedskin.yml`. Be sure that you are in your root hemoMIPs pipeline folder.

```bash
# create a new environment
conda env create -f envs/shedskin.yml -n shedskin

mkdir -p ~/miniconda3/envs/shedskin/etc/conda/activate.d
mkdir -p ~/miniconda3/envs/shedskin/etc/conda/deactivate.d

echo '#!/bin/sh
export LD_LIBRARY_PATH="$HOME/miniconda3/envs/shedskin/lib:$LD_LIBRARY_PATH"' > ~/miniconda3/envs/shedskin/etc/conda/activate.d/env_vars.sh

echo '#!/bin/sh
unset LD_LIBRARY_PATH' > ~/miniconda3/envs/shedskin/etc/conda/deactivate.d/env_vars.sh

source activate shedskin
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
The result should look like
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

Now we need to compile the `MergeTrimReads.py` script using Shed Skin:

```bash
# Go to the script folder
cd ../../scripts/pipeline2.0
# create Makefile and edit it
shedskin -e -L ~/miniconda3/envs/shedskin/include MergeTrimReads
sed -i '3s|$| -L ~/miniconda3/envs/shedskin/lib|' Makefile
# Compile!
make
cd ../../
```

### ensembl-vep

Finally, to run the pipeline Ensembl-VEP must be installed. In theory this could be done using conda environments as well but we found the vep conda recipe more laborious than the ensembl-vep installation guide:

https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html

Don't forget: to run ensembl-vep you also need to install the according homo-sapiens database, with the preinstalled ensembl-vep INSTALL.pl scripts.

We tested and ran the pipeline using ensembl-vep version 77 (https://github.com/Ensembl/ensembl-tools/archive/release/77.zip  )

### GATK v3

GATK v4 is included as a conda environment which automatically installs GATK v4.0.4.0 and all its dependencies.
If you prefer to run the pipeline using GATK v3 you need to install both GATK3 versions, if you trust GATK4, this step is not needed:

GATK 3.2.2 (https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.2-2-gec30cee)
GATK 3.4-46 (https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.4-46-gbc02625)


## Config

Almost ready to go:
We are aligning against the 1000 Genomes phase 2 build of the human reference: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/
You also need the bwa index of this file.
VEP uses following reference genome file: ftp://ftp.ensembl.org/pub/release-72/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.72.dna.toplevel.fa.gz
Also we add the CADD annotation cadd_v1.3 phase1_v3.20101123.vcf.gz to be found https://cadd.gs.washington.edu/download

Finally: Set the Paths and parameters in your config. You can find a template named "example_config.yml".

## Input

Put your NGS fastq files in input/dataset* together with:
- a MIP design file as generated by https://github.com/shendurelab/MIPGEN
- a sam header with the correct amount of samples
- a coresponding barcode sample assignment file
- the target coordinates of your MIPs
- the target coordinates of the captured sequences

example files can be found in the input/example_dataset folder

## Run pipeline

```bash
source activate hemoMIPs
# dry run to see if everything works
snakemake  --use-conda --configfile config.yml -n
# run the pipeline
snakemake  --use-conda --configfile config.yml
```
