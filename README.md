## _gvanno_ - *g*ermline *v*ariant *anno*tator

### Overview

The germline variant annotator (*gvanno*) is a simple, Docker-based software package intended for analysis and interpretation of human DNA variants of germline origin. Variants and genes are annotated with disease-related and functional associations from a wide range of sources (see below).

*gvanno* accepts query files encoded in the VCF format, and can analyze both SNVs and short InDels. The workflow relies heavily upon [Ensembl’s Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html), and [vcfanno](https://github.com/brentp/vcfanno). It produces an annotated VCF file and a file of tab-separated values (.tsv), the latter listing all annotations pr. variant record.

#### Annotation resources included in _gvanno_ - 1.1.0

* [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor v98 (GENCODE v31/v19 as the gene reference dataset)
* [dBNSFP](https://sites.google.com/site/jpopgen/dbNSFP) - Database of non-synonymous functional predictions (v4.0, May 2019)
* [gnomAD](http://gnomad.broadinstitute.org/) - Germline variant frequencies exome-wide (release 2.1, October 2018) - from VEP
* [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/) - Database of short genetic variants (build 152, January 2019) - from VEP
* [1000 Genomes Project - phase3](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) - Germline variant frequencies genome-wide (May 2013) - from VEP
* [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - Database of clinically related variants (November 2019)
* [DisGeNET](http://www.disgenet.org) - Database of gene-disease associations (v6.0, January 2019)
* [UniProt/SwissProt KnowledgeBase](http://www.uniprot.org) - Resource on protein sequence and functional information (2019_10, November 2019)
* [Pfam](http://pfam.xfam.org) - Database of protein families and domains (v32, Sept 2018)
* [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/home) - Catalog of published genome-wide association studies October 14th 2019)

### News
* November 22nd 2019 - **1.1.0 release**
     * Ability to install and run workflow using [Singularity](https://sylabs.io/docs/), excellent contribution by [@oskarvid](https://github.com/oskarvid), see step 1.1 in _Getting Started_
	* Data and software updates (ClinVar, UniProt, VEP)
* July 10th 2019 - **1.0.0 release**
     * Docker image update - VEP v97 (GENCODE 31/19)
     * Data bundle updates: ClinVar, UniProt, GWAS catalog
* May 21st 2019 - **0.9.0 release**
     * Data bundle updates: ClinVar, UniProt
	* Adding gene-disease associations from [Open Targets Platform](https://targetvalidation.org),([Carvalho-Silva et. al, NAR, 2019](https://www.ncbi.nlm.nih.gov/pubmed/30462303))
	* Moved *vcf-validation* configuration to command-line option

### Getting started

#### STEP 0: Python

An installation of Python (version _3.6_) is required to run *gvanno*. Check that Python is installed by typing `python --version` in your terminal window. In addition, a [Python library](https://github.com/uiri/toml) for parsing configuration files encoded with [TOML](https://github.com/toml-lang/toml) is needed. To install, simply run the following command:

   	pip install toml

#### STEP 1: Installation of Docker

1. [Install the Docker engine](https://docs.docker.com/engine/installation/) on your preferred platform
   - installing [Docker on Linux](https://docs.docker.com/engine/installation/linux/)
   - installing [Docker on Mac OS](https://docs.docker.com/engine/installation/mac/)
   - NOTE: We have not yet been able to perform enough testing on the Windows platform, and we have received feedback that particular versions of Docker/Windows do not work with PCGR (an example being [mounting of data volumes](https://github.com/docker/toolbox/issues/607))
2. Test that Docker is running, e.g. by typing `docker ps` or `docker images` in the terminal window
3. Adjust the computing resources dedicated to the Docker, i.e.:
   - Memory: minimum 5GB
   - CPUs: minimum 4
   - [How to - Mac OS X](https://docs.docker.com/docker-for-mac/#advanced)

##### STEP 1.1: Installation of Singularity (optional)
0. NB, this has only been tested with Singularity version 2.4.2, your mileage may vary with other versions.
1. [Install Singularity](https://sylabs.io/docs/)
2. Test that singularity works by running `singularity --version`
3. If you are in the gvanno directory, build the singularity image like so:

	`cd src`

  	Now you can build the singularity image with:

	`sudo ./buildSingularity.sh`

#### STEP 2: Download *gvanno* and data bundle

1. Download and unpack the [latest software release (1.1.0)](https://github.com/sigven/gvanno/releases/tag/v1.1.0)
2. Download and unpack the assembly-specific data bundle in the gvanno directory
   * [grch37 data bundle](https://drive.google.com/file/d/183S5XnTrPi1DDlVlO5_Cwmu73BAJZTAK) (approx 14Gb)
   * [grch38 data bundle](https://drive.google.com/file/d/1ZAjuPP7B2T0WgUsoZf20MTujZY9YNCN2) (approx 15Gb)
   * *Unpacking*: `gzip -dc gvanno.databundle.grch37.YYYYMMDD.tgz | tar xvf -`

    A _data/_ folder within the _gvanno-X.X_ software folder should now have been produced
3. Pull the [gvanno Docker image (1.1.0)](https://hub.docker.com/r/sigven/gvanno/) from DockerHub (approx 2Gb):
   * `docker pull sigven/gvanno:1.1.0` (gvanno annotation engine)

#### STEP 3: Input preprocessing

The *gvanno* workflow accepts a single input file:

  * An unannotated, single-sample VCF file (>= v4.2) with germline variants (SNVs/InDels)

We __strongly__ recommend that the input VCF is compressed and indexed using [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html). NOTE: If the input VCF contains multi-allelic sites, these will be subject to [decomposition](http://genome.sph.umich.edu/wiki/Vt#Decompose).

#### STEP 4: *gvanno* configuration

A few elements of the workflow can be figured using the *gvanno* configuration file (i.e. **gvanno.toml**), encoded in [TOML](https://github.com/toml-lang/toml) (an easy to read file format).

* Prediction of loss-of-function variants using VEP's LOFTEE plugin can be turned on in the configuration file (`lof_prediction = true`). Do note that this frequently increases the run time for VEP significantly.

#### STEP 5: Run example

Run the workflow with **gvanno.py**, which takes the following arguments and options:

		usage: gvanno.py [options] <QUERY_VCF> <GVANNO_DIR> <OUTPUT_DIR> <GENOME_ASSEMBLY>
		                           <CONFIG_FILE> <SAMPLE_ID> --container <docker|singularity>

		Germline variant annotation (gvanno) workflow for clinical and functional interpretation
		of germline nucleotide variants

		positional arguments:
		query_vcf           VCF input file with germline query variants (SNVs/InDels)
		gvanno_dir          gvanno base directory with accompanying data directory, e.g. ~/gvanno-1.0.0
		output_dir          Output directory
		{grch37,grch38}     grch37 or grch38
		configuration_file  gvanno configuration file (TOML format)
		sample_id           Sample identifier - prefix for output files
		--container         Run gvanno with docker or singularity

		optional arguments:
		-h, --help          show this help message and exit
		--force_overwrite   The script will fail with an error if the output file already exists. Force the overwrite of existing result files by using this flag
		--version           show program's version number and exit
		--no_vcf_validate   Skip validation of input VCF with Ensembl's vcf-validator


The _examples_ folder contains an example VCF file. Analysis of the example VCF can be performed by the following command:

`python ~/gvanno-1.1.0/gvanno.py ~/gvanno-1.1.0/examples/example_grch37.vcf.gz --container docker`
` ~/gvanno-1.1.0 ~/gvanno-1.1.0/examples grch37 ~/gvanno-1.1.0/gvanno.toml example`

This command will run the Docker-based *gvanno* workflow and produce the following output files in the _examples_ folder:

  1. __example_gvanno_pass_grch37.vcf.gz (.tbi)__ - Bgzipped VCF file with rich set of functional/clinical annotations
  2. __example_gvanno_pass_grch37.tsv.gz__ - Compressed TSV file with rich set of functional/clinical annotations

Similar files are produced for all variants, not only variants with a *PASS* designation in the VCF FILTER column.

### Documentation

Documentation of the various variant and gene annotations should be interrogated from the header of the annotated VCF file. The column names of the tab-separated values (TSV) file will be identical to the INFO tags that are documented in the VCF file.

### Contact

sigven@ifi.uio.no
