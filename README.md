# Nextflow Pipeline for Preprocessing GWAS Summary Statistics Files and Generating Files for the DyHealthNetLight Platform

Data preprocessing for the DyHealthNet Light platform is performed through a Nextflow-based pipeline that enables seamless deployment across different computing environments and automatically ensures scalability for large collections of GWAS summary statistics. 

The pipeline provides comprehensive functionality for preparing GWAS summary data for integration into the DyHealthNet Light platform, including data harmonization, variant annotation using the Ensembl Variant Effect Predictor, and, optionally, the execution of gene-based association analyses with the state-of-the-art tool MAGMA.

# Installations

## Java and Nextflow

Java and Nextflow need to be available on your machine. Details on how to install can be found here: https://www.nextflow.io/docs/latest/install.html.

## Conda Environment with VEP

Currently, VEP is not included in the main conda environment (envs/dyhealthnetlight_nf_pipeline.yml). Hence, you need to create an own conda environment and adapt the path of your conda in the configs/study_specific.config:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n ensembl-vep ensembl-vep=110.0
```

Note: The VEP version in conda (here 110.0) must match the VEP cache version specified by ensemblvep_cache_version in configs/study_specific.config !

# Reference Files

For the execution of MAGMA, reference files need to be downloaded. The 1000 Genome project is commonly used as a reference and both the reference data files and the gene locations of Ensembl can be downloaded from FUMA (https://fuma.ctglab.nl/downloadPage). 

We use these reference files because they are based on IDs consisting of chr:pos:alleles, instead of using rsIDs which are dependent on the dbSNP version.

Please then adapt the location of the reference files in the nextflow.config.


# Nextflow Config File

Study-specific configurations are required to be stored in configs/study-specific.config.

Computational parameters specifying the memory usage, time limits, etc. for the SLURM execution need to be specified in configs/compute_slurm.config.

The conda environment including the vep installation needs to be specified in nextflow.config.

# Run Nextflow Pipeline

```bash
nextflow run main.nf -profile slurm
```

# Replace Symlinks by Copying Files

Go to the output directory and execute this command:

```bash
find . -type l -exec sh -c '
  target=$(readlink -f "$1")
  if [ -f "$target" ]; then
    echo "Replacing symlink with real file: $1"
    cp --remove-destination "$target" "$1"
  else
    echo "Skipping broken symlink: $1 -> $target"
  fi
' _ {} \;
```


