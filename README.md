## *gvanno* - generic workflow for functional and clinical annotation of human DNA variants

<a href="https://github.com/sigven/gvanno/"><img src="src/gvanno_logo.png" align="right" height="130" width="113"/></a>

### Contents

-   [Overview](#overview)
-   [News](#news)
-   [Annotation resources](#annotation-resources)
-   [Getting started](#getting-started)
-   [Contact](#contact)

### Overview

The generic variant annotator (*gvanno*) is a software package intended for simple analysis and interpretation of human DNA variants. Variants and genes are annotated with disease-related and functional associations. Technically, the workflow is built with the [Docker](https://www.docker.com) technology, and it can also be installed through the [Singularity](https://sylabs.io/docs/) framework.

*gvanno* accepts query files encoded in the VCF format, and can analyze both SNVs and short insertions or deletions (indels). The workflow relies heavily upon [Ensembl's Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html), and [vcfanno](https://github.com/brentp/vcfanno). It produces an annotated VCF file and a file of tab-separated values (.tsv), the latter listing all annotations pr. variant record. Note that if your input VCF contains data (genotypes) from multiple samples (i.e. a multisample VCF), the output TSV file will contain one line/record **per sample variant**.

### News 

-   April 27th 2023 - **1.6.0 release**

    -   Added option `--oncogenicity_annotation` - classifies variants according to oncogenicity (Horak et al., Genet Med, 2022)
    -   Data updates: ClinVar, GENCODE, GWAS catalog, CancerMine
    -   Excluded Open Targets Platform from annotation

-   September 26th 2022 - **1.5.1 release**

    -   Added option `--vep_coding_only` - only report variants that fall into coding regions of transcripts (VEP option `--coding_only`)

-   September 24th 2022 - **1.5.0 release**

    -   Data updates: ClinVar, GENCODE, GWAS catalog, CancerMine, Open Targets Platform
    -   Software updates: VEP 107
    -   Excluded UniProt KB from annotation tracks

-   December 21st 2021 - **1.4.4 release**

    -   Data updates: ClinVar, GWAS catalog, CancerMine, UniProt KB, Open Targets Platform
    -   Software updates: VEP (v105)

-   August 25th 2021 - **1.4.3 release**

    -   Data updates: ClinVar, GWAS catalog, CancerMine, UniProt, Open Targets Platform

### Annotation resources (v1.6.0)

-   [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor v109 (GENCODE v43/v19 as the gene reference dataset)
-   [dBNSFP](https://sites.google.com/site/jpopgen/dbNSFP) - Database of non-synonymous functional predictions (v4.2, March 2021)
-   [gnomAD](http://gnomad.broadinstitute.org/) - Germline variant frequencies exome-wide (release 2.1, October 2018) - from VEP
-   [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/) - Database of short genetic variants (build 154) - from VEP
-   [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - Database of variants related to human health/disease phenotypes (April 2023)
-   [CancerMine](http://bionlp.bcgsc.ca/cancermine/) - literature-mined database of drivers, oncogenes and tumor suppressors in cancer (version 50, March 2023)
-   [Mutation hotspots](cancerhotspots.org) - Database of mutation hotspots in cancer
-   [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/home) - Catalog of published genome-wide association studies (March 27th 2023)

### Getting started

#### STEP 0: Prerequisites

-   *Python* 
    
    An installation of Python (version \>=*3.6*) is required to run *gvanno*. Check that Python is installed by typing `python --version` in your terminal window.

-   *Other utilities* 

    The downloading and unpacking of the reference data requires that the user has `bgzip` + `tabix` installed. See [here](http://www.htslib.org/download/) for instructions

#### STEP 1: Installation of Docker

1.  [Install the Docker engine](https://docs.docker.com/engine/installation/) on your preferred platform
    -   installing [Docker on Linux](https://docs.docker.com/engine/installation/linux/)
    -   installing [Docker on Mac OS](https://docs.docker.com/engine/installation/mac/)
    -   NOTE: We have not yet been able to perform enough testing on the Windows platform, and we have received feedback that particular versions of Docker/Windows do not work with gvanno (an example being [mounting of data volumes](https://github.com/docker/toolbox/issues/607))
2.  Test that Docker is running, e.g. by typing `docker ps` or `docker images` in the terminal window
3.  Adjust the computing resources dedicated to the Docker, i.e.:
    -   Memory: minimum 5GB
    -   CPUs: minimum 4
    -   [How to - Mac OS X](https://docs.docker.com/docker-for-mac/#advanced)

##### 1.1: Installation of Singularity (optional - in dev)

0.  **Note: this works for Singularity version 3.0 and higher**.

1.  [Install Singularity](https://sylabs.io/docs/)

2.  Test that singularity works by running `singularity --version`

3.  If you are in the gvanno directory, build the singularity image like so:

    `cd src`

    `sudo ./buildSingularity.sh`

#### STEP 2: Download *gvanno* and data bundle

1.  [Download the latest version](https://github.com/sigven/gvanno/releases/tag/v1.6.0)

2.  Download the assembly-specific VEP cache, and gvanno-specific reference data using the `download_gvanno_refdata.py` script, i.e.:

    -   `python download_gvanno_refdata.py --download_dir <PATH_TO_DOWNLOAD_DIR> --genome_assembly grch38`

    **NOTE**: This can take a considerable amount of time depending on your local bandwidth (approx 30Gb pr. assembly-specific bundle)

3.  Pull the [gvanno Docker image (1.6.0)](https://hub.docker.com/r/sigven/gvanno/) from DockerHub (approx 4.11Gb):

    -   `docker pull sigven/gvanno:1.6.0` (gvanno annotation engine)

#### STEP 3: Input preprocessing

The *gvanno* workflow accepts a single input file:

-   An unannotated, single-sample VCF file (\>= v4.2) with germline variants (SNVs/InDels)

We **strongly** recommend that the input VCF is compressed and indexed using [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html). NOTE: If the input VCF contains multi-allelic sites, these will be subject to [decomposition](http://genome.sph.umich.edu/wiki/Vt#Decompose).

#### STEP 5: Run example

Run the workflow with **gvanno.py**, which takes the following arguments and options:

```         
usage:
gvanno.py -h [options]
--query_vcf <QUERY_VCF>
--gvanno_dir <GVANNO_DIR>
--output_dir <OUTPUT_DIR>
--genome_assembly <grch37|grch38>
--sample_id <SAMPLE_ID>
--container <docker|singularity>

gvanno - workflow for functional and clinical annotation of germline nucleotide variants

Required arguments:
--query_vcf QUERY_VCF
                VCF input file with germline query variants (SNVs/InDels).
--gvanno_dir GVANNO_DIR
                Directory that contains the gvanno data bundle, e.g. ~/gvanno-1.6.0
--output_dir OUTPUT_DIR
                Output directory
--genome_assembly {grch37,grch38}
                Genome assembly build: grch37 or grch38
--container {docker,singularity}
                Run gvanno with docker or singularity
--sample_id SAMPLE_ID
                Sample identifier - prefix for output files

VEP optional arguments:
--vep_regulatory      Enable Variant Effect Predictor (VEP) to look for overlap with regulatory regions (option --regulatory in VEP).
--vep_gencode_all     Consider all GENCODE transcripts with Variant Effect Predictor (VEP) (option --gencode_basic in VEP is used by default in gvanno).
--vep_lof_prediction  Predict loss-of-function variants with Loftee plugin in Variant Effect Predictor (VEP), default: False
--vep_n_forks VEP_N_FORKS
                Number of forks for Variant Effect Predictor (VEP) processing, default: 4
--vep_buffer_size VEP_BUFFER_SIZE
                Variant buffer size (variants read into memory simultaneously) for Variant Effect Predictor (VEP) processing
                - set lower to reduce memory usage, higher to increase speed, default: 500
--vep_pick_order VEP_PICK_ORDER
                Comma-separated string of ordered transcript properties for primary variant pick in
                Variant Effect Predictor (VEP) processing, default: canonical,appris,biotype,ccds,rank,tsl,length,mane
--vep_skip_intergenic
                Skip intergenic variants in Variant Effect Predictor (VEP) processing, default: False
--vep_coding_only
          Only report variants falling into coding regions of transcripts (VEP), default: False

Other optional arguments:
--force_overwrite     By default, the script will fail with an error if any output file already exists.
                You can force the overwrite of existing result files by using this flag, default: False
--version             show program's version number and exit
--no_vcf_validate     Skip validation of input VCF with Ensembl's vcf-validator, default: False
--docker_uid DOCKER_USER_ID
                Docker user ID. default is the host system user ID. If you are experiencing permission errors, try setting this up to root (`--docker-uid root`)
--vcfanno_n_processes VCFANNO_N_PROCESSES
                Number of processes for vcfanno processing (see https://github.com/brentp/vcfanno#-p), default: 4
--oncogenicity_annotation
                    Classify variants according to oncogenicity (Horak et al., Genet Med, 2022)
--debug             Print full Docker/Singularity commands to log and do not delete intermediate files with warnings etc.
```

The *examples* folder contains an example VCF file. Analysis of the example VCF can be performed by the following command:

```         
python ~/gvanno-1.6.0/gvanno.py
--query_vcf ~/gvanno-1.6.0/examples/example.grch37.vcf.gz
--gvanno_dir ~/gvanno-1.6.0
--output_dir ~/gvanno-1.6.0
--sample_id example
--genome_assembly grch37
--container docker
--force_overwrite
```

This command will run the Docker-based *gvanno* workflow and produce the following output files in the *examples* folder:

1.  **example_gvanno_pass_grch37.vcf.gz (.tbi)** - Bgzipped VCF file with rich set of functional/clinical annotations
2.  **example_gvanno_pass_grch37.tsv.gz** - Compressed TSV file with rich set of functional/clinical annotations

Similar files are produced for all variants, not only variants with a *PASS* designation in the VCF FILTER column.

### Documentation

Documentation of the various variant and gene annotations should be interrogated from the header of the annotated VCF file. The column names of the tab-separated values (TSV) file will be identical to the INFO tags that are documented in the VCF file.

### Contact

sigven AT ifi.uio.no
