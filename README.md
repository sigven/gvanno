## _gvanno_ - workflow for functional and clinical annotation of germline nucleotide variants

### Contents

- [Overview](#overview)
- [News](#news)
- [Annotation resources](#annotation-resources)
- [Getting started](#getting-started)
- [Contact](#contact)

### Overview

The germline variant annotator (*gvanno*) is a simple software package intended for analysis and interpretation of human DNA variants of germline origin. Variants and genes are annotated with disease-related and functional associations from a wide range of sources (see below). Technically, the workflow is built with the [Docker](https://www.docker.com) technology, and it can also be installed through the [Singularity](https://sylabs.io/docs/) framework.

*gvanno* accepts query files encoded in the VCF format, and can analyze both SNVs and short InDels. The workflow relies heavily upon [Ensembl’s Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html), and [vcfanno](https://github.com/brentp/vcfanno). It produces an annotated VCF file and a file of tab-separated values (.tsv), the latter listing all annotations pr. variant record. Note that if your input VCF contains data (genotypes) from multiple samples (i.e. a multisample VCF), the output TSV file will contain one line/record __per sample variant__.

### News
* December 7th 2020 - **1.4.1 release**
  * Data updates (ClinVar, UniProt, GWAS Catalog, Open Targets Platform)
  * Software update (VEP 102)
  * Skipped DisGenet annotations (Open Targets serve similar purpose)
* September 29th 2020 - **1.4.0 release**
  * Data updates (ClinVar, UniProt, GWAS Catalog, Open Targets Platform)
  * Software updates (VEP 101)
  * Configuration through TOML file is omitted - all configurations are now encoded as optional arguments to the main Python script (`gvanno.py`)
* June 30th 2020 - **1.3.2 release**
     * Data updates (ClinVar, UniProt, GWAS Catalog, Open Targets Platform, Pfam, dbNSFP)
	     * Using GENCODE v34 as the correct transcript assembly for grch38 (see [issue](https://github.com/Ensembl/ensembl-vep/issues/749))
		* Three new variant effect predictions from dbNSFP added: [ClinPred](https://doi.org/10.1016/j.ajhg.2018.08.005), [LIST-S2](https://doi.org/10.1093/nar/gkaa288), and [BayesDel](https://doi.org/10.1002/humu.23158)
	* Added VEP plugin [NearestExonJB](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#plugins_existing)
	     * Annotates relative position (to the exon-intron junction) of variants in introns and exons (fields in output: INTRON_POSITION, EXON_POSITION)
* May 8th 2020 - **1.3.0 release**
     * Upgrade of VEP (v100) - GENCODE release 33 (grch38)
	* Data updates (ClinVar, UniProt, GWAS Catalog, Open Targets Platform)
* November 22nd 2019 - **1.1.0 release**
     * Ability to install and run workflow using [Singularity](https://sylabs.io/docs/), excellent contribution by [@oskarvid](https://github.com/oskarvid), see step 1.1 in _Getting Started_
	* Data and software updates (ClinVar, UniProt, VEP)


### Annotation resources

* [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor v102 (GENCODE v36/v19 as the gene reference dataset)
* [dBNSFP](https://sites.google.com/site/jpopgen/dbNSFP) - Database of non-synonymous functional predictions (v4.1, June 2020)
* [gnomAD](http://gnomad.broadinstitute.org/) - Germline variant frequencies exome-wide (release 2.1, October 2018) - from VEP
* [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/) - Database of short genetic variants (build 153) - from VEP
* [1000 Genomes Project - phase3](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) - Germline variant frequencies genome-wide (May 2013) - from VEP
* [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - Database of clinically related variants (December 2020)
* [Open Targets Platform](https://targetvalidation.org) - Target-disease and target-drug associations (2020_11, November 2020)
* [UniProt/SwissProt KnowledgeBase](http://www.uniprot.org) - Resource on protein sequence and functional information (2020_06, December 2020)
* [Pfam](http://pfam.xfam.org) - Database of protein families and domains (v33.1, May 2020)
* [NHGRI-EBI GWAS Catalog](https://www.ebi.ac.uk/gwas/home) - Catalog of published genome-wide association studies (December 2nd 2020)


### Getting started

#### STEP 0: Python

An installation of Python (version _3.6_) is required to run *gvanno*. Check that Python is installed by typing `python --version` in your terminal window.

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

##### 1.1: Installation of Singularity (optional)

0. **Note: this has only been tested with Singularity version 2.4.2, your mileage may vary with other versions**.
1. [Install Singularity](https://sylabs.io/docs/)
2. Test that singularity works by running `singularity --version`
3. If you are in the gvanno directory, build the singularity image like so:

	`cd src`

	`sudo ./buildSingularity.sh`

#### STEP 2: Download *gvanno* and data bundle

1. Download and unpack the [latest software release (1.4.1)](https://github.com/sigven/gvanno/releases/tag/v1.4.1)
2. Download and unpack the assembly-specific data bundle in the gvanno directory
   * [grch37 data bundle](http://insilico.hpc.uio.no/pcgr/gvanno/gvanno.databundle.grch37.20201206.tgz) (approx 16Gb)
   * [grch38 data bundle](http://insilico.hpc.uio.no/pcgr/gvanno/gvanno.databundle.grch38.20201206.tgz) (approx 17Gb)
   * *Unpacking*: `gzip -dc gvanno.databundle.grch37.YYYYMMDD.tgz | tar xvf -`

    A _data/_ folder within the _gvanno-X.X_ software folder should now have been produced
3. Pull the [gvanno Docker image (1.4.1)](https://hub.docker.com/r/sigven/gvanno/) from DockerHub (approx 2.3Gb):
   * `docker pull sigven/gvanno:1.4.1` (gvanno annotation engine)

#### STEP 3: Input preprocessing

The *gvanno* workflow accepts a single input file:

  * An unannotated, single-sample VCF file (>= v4.2) with germline variants (SNVs/InDels)

We __strongly__ recommend that the input VCF is compressed and indexed using [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html). NOTE: If the input VCF contains multi-allelic sites, these will be subject to [decomposition](http://genome.sph.umich.edu/wiki/Vt#Decompose).

#### STEP 5: Run example

Run the workflow with **gvanno.py**, which takes the following arguments and options:

	usage:
	gvanno.py -h [options]
	--query_vcf QUERY_VCF
	--gvanno_dir GVANNO_DIR
	--output_dir OUTPUT_DIR
	--genome_assembly grch37|grch38
	--sample_id SAMPLE_ID
	--container docker|singularity

	gvanno - workflow for functional and clinical annotation of germline nucleotide variants

	Required arguments:
	--query_vcf QUERY_VCF
			    VCF input file with germline query variants (SNVs/InDels).
	--gvanno_dir GVANNO_DIR
			    Directory that contains the gvanno data bundle, e.g. ~/gvanno-1.4.1
	--output_dir OUTPUT_DIR
			    Output directory
	--genome_assembly {grch37,grch38}
			    Genome assembly build: grch37 or grch38
	--container {docker,singularity}
			    Run gvanno with docker or singularity
	--sample_id SAMPLE_ID
			    Sample identifier - prefix for output files

	Optional arguments:
	--force_overwrite     By default, the script will fail with an error if any output file already exists.
			    You can force the overwrite of existing result files by using this flag, default: False
	--version             show program's version number and exit
	--no_vcf_validate     Skip validation of input VCF with Ensembl's vcf-validator, default: False
	--lof_prediction      Predict loss-of-function variants with Loftee plugin in Variant Effect Predictor (VEP), default: False
	--vep_n_forks VEP_N_FORKS
			    Number of forks for Variant Effect Predictor (VEP) processing, default: 4
	--vep_buffer_size VEP_BUFFER_SIZE
			    Variant buffer size (variants read into memory simultaneously) for Variant Effect Predictor (VEP) processing
			    - set lower to reduce memory usage, default: 5000
	--vep_pick_order VEP_PICK_ORDER
			    Comma-separated string of ordered transcript properties for primary variant pick in
			    Variant Effect Predictor (VEP) processing, default: canonical,appris,biotype,ccds,rank,tsl,length,mane
	--vep_skip_intergenic
			    Skip intergenic variants in Variant Effect Predictor (VEP) processing, default: False
	--vcfanno_n_processes VCFANNO_N_PROCESSES
			    Number of processes for vcfanno processing (see https://github.com/brentp/vcfanno#-p), default: 4


The _examples_ folder contains an example VCF file. Analysis of the example VCF can be performed by the following command:

	python ~/gvanno-1.4.1/gvanno.py
	--query_vcf ~/gvanno-1.4.1/examples/example.grch37.vcf.gz
	--gvanno_dir ~/gvanno-1.4.1
	--output_dir ~/gvanno-1.4.1
	--sample_id example
	--genome_assembly grch37
	--container docker
	--force_overwrite

This command will run the Docker-based *gvanno* workflow and produce the following output files in the _examples_ folder:

  1. __example_gvanno_pass_grch37.vcf.gz (.tbi)__ - Bgzipped VCF file with rich set of functional/clinical annotations
  2. __example_gvanno_pass_grch37.tsv.gz__ - Compressed TSV file with rich set of functional/clinical annotations

Similar files are produced for all variants, not only variants with a *PASS* designation in the VCF FILTER column.

### Documentation

Documentation of the various variant and gene annotations should be interrogated from the header of the annotated VCF file. The column names of the tab-separated values (TSV) file will be identical to the INFO tags that are documented in the VCF file.

### Contact

sigven AT ifi.uio.no
