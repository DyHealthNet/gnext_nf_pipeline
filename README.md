# Nextflow Preprocessing Pipeline for DyHealthNet Light Platform

Data preprocessing for the DyHealthNet Light platform is performed through a Nextflow-based pipeline that enables seamless deployment across different computing environments and automatically ensures scalability for large collections of GWAS summary statistics. 

The pipeline provides comprehensive functionality for preparing GWAS summary data for integration into the DyHealthNet Light platform, including data harmonization, variant annotation using the Ensembl Variant Effect Predictor, and, optionally, the execution of gene-based association analyses with the state-of-the-art tool MAGMA.

<img width="3308" height="1480" alt="Nextflow_pipeline" src="https://github.com/user-attachments/assets/b607dcf8-7afe-48d9-a306-64fd37bb97ba" />

# Installations

Execution of the pipeline requires the installation of Java, Nextflow, and Conda.

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

## Study-Specific Configs

| **Parameter** | **Example Value** | **Description** |
|----------------|-----------------------------|-----------------|
| `pheno_file` | `/storage03/larend/lipids_data/lipids_case_study_traits.csv` | Path to the phenotype configuration file containing trait names and corresponding GWAS summary statistic file paths. |
| `out_dir` | `/storage03/larend/lipids_data/lipids_nf_out` | Output directory where all processed results and intermediate files will be stored. |
| `pheno_batch_size` | `5` | Number of phenotypes processed in parallel within a single batch. |
| `steps` | `["gene_statistics", "gwas_exploration"]` | List of workflow steps to be executed (e.g., gene-based statistics, exploratory analysis). |
| `chr_column` | `3` | Column index (1-based) of the chromosome field in the GWAS summary statistics file. |
| `pos_column` | `4` | Column index of the variant position field. |
| `ref_column` | `5` | Column index of the reference allele field. |
| `alt_column` | `6` | Column index of the alternate allele field. |
| `pval_column` | `16` | Column index of the p-value field. |
| `beta_column` | `9` | Column index of the effect size (beta) field. |
| `se_column` | `10` | Column index of the standard error field. |
| `af_column` | `8` | Column index of the allele frequency field. |
| `pval_neglog10` | `false` | Indicates whether p-values are stored as negative log10 values (`true`) or raw p-values (`false`). |
| `ensemblvep_species` | `'homo_sapiens'` | Species identifier for Ensembl VEP annotation. |
| `ensemblvep_genome` | `'GRCh37'` | Genome assembly version used for annotation. |
| `ensemblvep_cache_version` | `110` | Version of the Ensembl VEP cache used during annotation. Must match the installed VEP version. |
| `ensemblvep_conda_env` | `/home/larend/miniforge3/envs/vep` | Path to the Conda environment containing the VEP installation. |
| `magma_reference_plink` | `/storage03/larend/MAGMA_ref/FUMA/EUR/EUR` | Path to the PLINK reference panel used by MAGMA for gene-based analyses. |
| `magma_window_up` | `10` | Upstream window size (in kb) used when mapping variants to genes in MAGMA. |
| `magma_window_down` | `10` | Downstream window size (in kb) used when mapping variants to genes in MAGMA. |
| `magma_gene_location` | `/storage03/larend/MAGMA_ref/FUMA/ENSGv110.coding.genes.txt` | Path to the Ensembl gene location file used for MAGMA analyses. |

## Slurm Configs

| **Parameter** | **Example / Default Value** | **Description** |
|----------------|-----------------------------|-----------------|
| `global_maxForks` | `10` | Defines the global maximum number of processes that can run in parallel across the workflow. |
| `normalize_cpus` | `12` | Number of CPU cores allocated for the normalization step. |
| `normalize_memory` | `'64GB'` | Memory allocated for the normalization step. |
| `vcf_cpus` | `16` | Number of CPU cores allocated for variant annotation and processing (e.g., VEP execution). |
| `vcf_memory` | `'64GB'` | Memory allocated for VCF processing and annotation. |
| `manhattan_qq_cpus` | `12` | Number of CPU cores used for generating Manhattan and QQ plots. |
| `manhattan_qq_memory` | `'64GB'` | Memory allocated for the Manhattan and QQ plot generation step. |
| `chrom_bgz_maxForks` | `23` | Maximum number of chromosome compression tasks (`bgzip`) that can run simultaneously. |
| `chrom_bgz_memory` | `'64GB'` | Memory allocated for per-chromosome compression tasks. |
| `chrom_bgz_cpus` | `8` | Number of CPU cores allocated for per-chromosome compression (`bgzip`) tasks. |
| `magma_cpus` | `32` | Number of CPU cores allocated for MAGMA gene-based analyses. |
| `magma_memory` | `'64GB'` | Memory allocated for MAGMA analyses. |
| `magma_input_cpus` | `32` | Number of CPU cores allocated for preparing MAGMA input files. |
| `magma_input_memory` | `'64GB'` | Memory allocated for MAGMA input file preparation. |

| **Setting** | **Value** | **Description** |
|--------------|------------|-----------------|
| `process.conda` | `/home/larend/miniforge3/envs/dyhealthnetlight_nf_pipeline` | Path to the Conda environment used for all SLURM-executed processes. |
| `process.executor` | `"slurm"` | Defines the workflow executor. Tasks are submitted to a SLURM cluster. |
| `process.queue` | `"slow-mc2"` | Specifies the SLURM queue (partition) to which jobs are submitted. |
| `process.maxForks` | `params.global_maxForks` (default `8`) | Maximum number of concurrent processes allowed under the SLURM profile. |

## Base Config

These parameters are typically not intended to be modified.

| **Parameter** | **Example / Default Value** | **Description** |
|----------------|-----------------------------|-----------------|
| `manhattan_num_unbinned` | `500` | Number of unbinned variants displayed in the Manhattan plot to preserve the most significant points without binning. |
| `manhattan_peak_max_count` | `500` | Maximum number of peaks (significant loci) displayed in the Manhattan plot for readability and performance. |
| `manhattan_peak_pval_threshold` | `1e-6` | P-value significance threshold used for identifying peaks in the Manhattan plot. Variants below this value are considered significant. |
| `manhattan_peak_sprawl_dist` | `200_000` | Minimum genomic distance (in base pairs) required to distinguish separate peaks in the Manhattan plot. Peaks closer than this are merged. |
| `top_hits_pval_cutoff` | `1e-6` | P-value threshold used to select top associated variants for downstream analysis. |
| `top_hits_max_limit` | `10,000` | Maximum number of top associated variants reported or exported after filtering by p-value. |
| `ensemblvep_distance_up` | `5,000` | Upstream distance (in base pairs) used by Ensembl VEP when mapping variants to nearby genes. |
| `ensemblvep_distance_down` | `5,000` | Downstream distance (in base pairs) used by Ensembl VEP when mapping variants to nearby genes. |


# Run Nextflow Pipeline

Once all parameters have been correctly specified, the pipeline can be executed.

```bash
nextflow run main.nf -profile slurm
```

# Replace Symlinks by Copying Files

The Nextflow pipeline generates symbolic links in the output directory. Therefore, as a final step, navigate to the output directory and execute the following command:

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


