## _gvanno_ - *g*ermline *v*ariant *anno*tator

### Overview

The germline variant annotator (*gvanno*) is a simple, stand-alone software package intended for analysis and interpretation of human germline calls. It accepts query files encoded in the VCF format, and can analyze both SNVs and short InDels. The software extends basic annotations from [Ensembl’s Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html) with up-to-date functional and clinical variant annotations retrieved flexibly through [vcfanno](https://github.com/brentp/vcfanno). The workflow produces an annotated VCF and a file of tab-separated values (.tsv), the latter listing all annotations pr. variant record.

#### Annotation resources included in _gvanno_ - 0.2.0

* [VEP v90](http://www.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor release 90 (GENCODE v27 as the gene reference dataset)
* [dBNSFP v3.4](https://sites.google.com/site/jpopgen/dbNSFP) - Database of non-synonymous functional predictions (March 2017)
* [gnomAD r1](http://gnomad.broadinstitute.org/) - Germline variant frequencies exome-wide (March 2017)
* [dbSNP b147](http://www.ncbi.nlm.nih.gov/SNP/) - Database of short genetic variants (April 2016)
* [1000 Genomes Project - phase3](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) - Germline variant frequencies genome-wide (May 2013)
* [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - Database of clinically related variants (November 2017)
* [DoCM](http://docm.genome.wustl.edu) - Database of curated mutations (v3.2, April 2016)
* [CIViC](http://civic.genome.wustl.edu) - Clinical interpretations of variants in cancer (November 11th 2017)
* [DisGeNET](http://www.disgenet.org) - Database of gene-disease associations (May 2017)
* [UniProt/SwissProt KnowledgeBase 2017_10](http://www.uniprot.org) - Resource on protein sequence and functional information (October 2017)
* [Pfam v31](http://pfam.xfam.org) - Database of protein families and domains (March 2017)
* [TSGene v2.0](http://bioinfo.mc.vanderbilt.edu/TSGene/) - Tumor suppressor/oncogene database (November 2015)

### Documentation

Coming

### Getting started

#### STEP 0: Python

A local installation of Python (it has been tested with [version 2.7.13](https://www.python.org/downloads/)) is required to run gvanno. Check that Python is installed by typing `python --version` in a terminal window. In addition, a [Python library](https://github.com/uiri/toml) for parsing configuration files encoded with [TOML](https://github.com/toml-lang/toml) is needed. To install, simply run the following command:

   	pip install toml

#### STEP 1: Installation of Docker

1. [Install the Docker engine](https://docs.docker.com/engine/installation/) on your preferred platform
   - installing [Docker on Linux](https://docs.docker.com/engine/installation/linux/)
   - installing [Docker on Mac OS](https://docs.docker.com/engine/installation/mac/)
   - NOTE: We have not yet been able to perform enough testing on the Windows platform, and we have received feedback that particular versions of Docker/Windows do not work with _gvanno_ (an example being [mounting of data volumes](https://github.com/docker/toolbox/issues/607))
2. Test that Docker is running, e.g. by typing `docker ps` or `docker images` in the terminal window
3. Adjust the computing resources dedicated to the Docker, i.e.:
   - Memory: minimum 5GB
   - CPUs: minimum 4
   - [How to - Mac OS X](https://docs.docker.com/docker-for-mac/#advanced)

#### STEP 2: Download _gvanno_

1. Download and unpack the [latest software release (0.2.0)](https://github.com/sigven/gvanno/releases/tag/v0.2.0)
2. Download and unpack the data bundle (approx. 15Gb) in the _gvanno_ directory
   * Download [the accompanying data bundle](https://drive.google.com/file/d/1uDFanR2LURgDjO_EB0ADzWBp5myE2rqn/) from Google Drive to `~/gvanno-X.X` (replace _X.X_ with the version number, e.g `~/gvanno-0.2.0`)
   * Unpack the data bundle, e.g. through the following Unix command: `gzip -dc gvanno.databundle.GRCh37.YYYYMMDD.tgz | tar xvf -`
    A _data/_ folder within the _gvanno-X.X_ software folder should now have been produced
3. Pull the [_gvanno_ Docker image (0.2.0)](https://hub.docker.com/r/sigven/gvanno/) from DockerHub (approx 4.2Gb):
   * `docker pull sigven/gvanno:0.2.0` (_gvanno_ annotation engine)

#### STEP 3: Input preprocessing

The _gvanno_ workflow accepts a single input file:

  * An unannotated (preferably single sample) VCF file (>= v4.2) with called germline variants (SNVs/InDels)
  * __NOTE__: GRCh37 is currently supported as the reference genome build
  * We __strongly__ recommend that the input VCF is compressed and indexed using [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html)
  * If the input VCF contains multi-allelic sites, these will be subject to [decomposition](http://genome.sph.umich.edu/wiki/Vt#Decompose)


#### STEP 4: Run example

Run the workflow with **gvanno.py**, which takes the following arguments and options:

	usage: gvanno.py [-h] [--input_vcf INPUT_VCF]
			 [--force_overwrite] [--version]
			 gvanno_dir output_dir configuration_file sample_id

	Germline variant annotation workflow for clinical and functional interpretation of
	single nucleotide variants and short insertions/deletions

	positional arguments:
	gvanno_dir              gvanno base directory with accompanying data directory,
					e.g. ~/gvanno-0.2.0
	output_dir            Output directory
	configuration_file    gvanno configuration file (TOML format)
	sample_id             Sample identifier - prefix for
					output files

	optional arguments:
	-h, --help            show this help message and exit
	--input_vcf INPUT_VCF
					VCF input file with somatic query variants
					(SNVs/InDels). Note: GRCh37 is currently the only
					reference genome build supported (default: None)
	--force_overwrite     By default, the script will fail with an error if any
					output file already exists. You can force the
					overwrite of existing result files by using this flag
					(default: False)
	--version             show program's version number and exit




### Contact

sigven@ifi.uio.no
