## _gvanno_ - *g*ermline *v*ariant *anno*tator

### Overview

The germline variant annotator (*gvanno*) is a simple, Docker-based software package intended for analysis and interpretation of human DNA variants of germline origin. It accepts query files encoded in the VCF format, and can analyze both SNVs and short InDels. The workflow is largely based on [Ensembl’s Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html), and [vcfanno](https://github.com/brentp/vcfanno). It produces an annotated VCF file and a file of tab-separated values (.tsv), the latter listing all annotations pr. variant record.

#### Annotation resources included in _gvanno_ - 0.3.1

* [VEP v92](http://www.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor release 92 (GENCODE v19/v28 as the gene reference dataset)
* [dBNSFP v3.5](https://sites.google.com/site/jpopgen/dbNSFP) - Database of non-synonymous functional predictions (August 2017)
* [gnomAD r2](http://gnomad.broadinstitute.org/) - Germline variant frequencies exome-wide (February 2017) - from VEP
* [dbSNP b150](http://www.ncbi.nlm.nih.gov/SNP/) - Database of short genetic variants (February 2017) - from VEP
* [1000 Genomes Project - phase3](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) - Germline variant frequencies genome-wide (May 2013) - from VEP
* [ClinVar 20180603](http://www.ncbi.nlm.nih.gov/clinvar/) - Database of clinically related variants (June 2018)
* [DisGeNET](http://www.disgenet.org) - Database of gene-disease associations (v5.0, May 2017)
* [UniProt/SwissProt KnowledgeBase 2018_06](http://www.uniprot.org) - Resource on protein sequence and functional information (June 2018)
* [Pfam v31](http://pfam.xfam.org) - Database of protein families and domains (March 2017)
* [TSGene v2.0](http://bioinfo.mc.vanderbilt.edu/TSGene/) - Tumor suppressor/oncogene database (November 2015)

### News


* July 5th 2018 - **0.3.1 release**
     * Data bundle updates (ClinVar, UniProt)
     * Addition of [VEP LofTee plugin](https://github.com/konradjk/loftee) - predicts loss-of-function variants
* April 20th 2018 - **0.3.0 release**
	* Runs under Python3
	* VEP version 92
	* Support for grch38
	* Data bundle updates (ClinVar, UniProt)

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

#### STEP 2: Download *gvanno* and data bundle

1. Download and unpack the [latest software release (0.3.1)](https://github.com/sigven/gvanno/releases/tag/v0.3.1)
2. Download and unpack the assembly-specific data bundle in the PCGR directory
   * [grch37 data bundle](https://drive.google.com/file/d/15NbYwwnb8J5IGhL6-RJXpAeQ-xqzjc5F/) (approx 9Gb)
   * [grch38 data bundle](https://drive.google.com/file/d/1hr4MShsEh2Xf-_bBgDPi7t-vj32XrWJ0/) (approx 9Gb)
   * *Unpacking*: `gzip -dc gvanno.databundle.grch37.YYYYMMDD.tgz | tar xvf -`

    A _data/_ folder within the _gvanno-X.X_ software folder should now have been produced
3. Pull the [gvanno Docker image (0.3.1)](https://hub.docker.com/r/sigven/gvanno/) from DockerHub (approx 2.5Gb):
   * `docker pull sigven/gvanno:0.3.1` (gvanno annotation engine)

#### STEP 3: Input preprocessing

The *gvanno* workflow accepts a single input file:

  * An unannotated, single-sample VCF file (>= v4.2) with germline variants (SNVs/InDels)

We __strongly__ recommend that the input VCF is compressed and indexed using [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html). NOTE: If the input VCF contains multi-allelic sites, these will be subject to [decomposition](http://genome.sph.umich.edu/wiki/Vt#Decompose).

#### STEP 4: *gvanno* configuration

A few elements of the workflow can be figured using the *gvanno* configuration file, encoded in [TOML](https://github.com/toml-lang/toml) (an easy to read file format).

The initial step of the workflow performs [VCF validation](https://github.com/EBIvariation/vcf-validator) on the input VCF file. This procedure is very strict, and often causes the workflow to return an error due to various violations of the VCF specification. If the user trusts that the most critical parts of the input VCF is properly encoded,  a setting in the configuration file (`vcf_validation = false`) can be used to turn off VCF validation.

#### STEP 5: Run example

Run the workflow with **gvanno.py**, which takes the following arguments and options:

	usage: gvanno.py [-h] [--input_vcf INPUT_VCF] [--force_overwrite] [--version]
			  gvanno_dir output_dir {grch37,grch38} configuration_file
			  sample_id

	Germline variant annotation (gvanno) workflow for clinical and functional
	interpretation of germline nucleotide variants

	positional arguments:
	gvanno_dir            gvanno base directory with accompanying data
				    directory, e.g. ~/gvanno-0.3.1
	output_dir            Output directory
	{grch37,grch38}       grch37 or grch38
	configuration_file    gvanno configuration file (TOML format)
	sample_id             Sample identifier - prefix for output files

	optional arguments:
	-h, --help            show this help message and exit
	--input_vcf INPUT_VCF
				    VCF input file with somatic query variants
				    (SNVs/InDels) (default: None)
	--force_overwrite     The script will fail with an error if the output file
				    already exists. Force the overwrite of existing result
				    files by using this flag (default: False)
	--version             show program's version number and exit


The _examples_ folder contains an example VCF file. It also contains a *gvanno* configuration file. Analysis of the example VCF can be performed by the following command:

`python ~/gvanno-0.3.1/gvanno.py --input_vcf ~/gvanno-0.3.1/examples/example.vcf.gz`
` ~/gvanno-0.3.1 ~/gvanno-0.3.1/examples grch37 ~/gvanno-0.3.1/examples/gvanno_config.toml example`


This command will run the Docker-based *gvanno* workflow and produce the following output files in the _examples_ folder:

  1. __example_gvanno_pass_grch37.vcf.gz (.tbi)__ - Bgzipped VCF file with rich set of functional/clinical annotations
  2. __example_gvanno_pass_grch37.tsv.gz__ - Compressed TSV file with rich set of functional/clinical annotations

Similar files are produced for all variants, not only variants with a *PASS* designation in the VCF FILTER column.

Documentation of the various variant and gene annotations should be interrogated from the header of the annotated VCF file.



### Contact

sigven@ifi.uio.no
