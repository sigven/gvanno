## gvanno - germline variant annotator

### Overview

The germline variant annotator (gvanno) is a stand-alone software package intended for analysis and interpretation of human germline calls. It accepts query files encoded in the VCF format, and can analyze both SNVs and short InDels. The software extends basic gene and variant annotations from the [Ensembl’s Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html) with up-to-date annotations retrieved flexibly through [vcfanno](https://github.com/brentp/vcfanno) (Pedersen et al., Genome Biology, 2016).

#### Annotation resources included in gvanno - v0.1

* [VEP v85](http://www.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor release 85 (GENCODE v19 as the gene model)
* [dBNSFP v3.2](https://sites.google.com/site/jpopgen/dbNSFP) - Database of non-synonymous functional predictions (March 2016)
* [ExAC r0.3.1](http://exac.broadinstitute.org/) - Germline variant frequencies exome-wide (March 2016)
* [dbSNP b147](http://www.ncbi.nlm.nih.gov/SNP/) - Database of short genetic variants (April 2016)
* [1000Genomes phase3](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) - Germline variant frequencies genome-wide (May 2013)
* [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - Database of clinically related variants (Nov 2016)
* [DoCM](http://docm.genome.wustl.edu) - Database of curated mutations (v3.2, April 2016)
* [UniProt/SwissProt KnowledgeBase 2016_09](http://www.uniprot.org) - Resource on protein sequence and functional information (Sep 2016)
* [Pfam v30](http://pfam.xfam.org) - Database of protein families and domains (June 2016)
* [TSGene v2.0](http://bioinfo.mc.vanderbilt.edu/TSGene/) - Tumor suppressor/oncogene database (November 2015)
* [DisGenNet v4.0 - gene-disease associations](http://www.disgenet.org) (April 2016)

### Getting started

#### STEP 1: Installation of Docker

1. TODO (Ghis): Bullet-proof Docker installation instructions (Mac, Windows(?), Linux)
2. TODO (Ghis): Test that Docker is running
2. TODO (Ghis): Recommended settings for Docker VM (Memory: 5GB, CPUs: 4), how to set them

#### STEP 2: Installation of gvanno (GRCh37)

1. Make a gvanno directory
2. Download and unpack the data bundle (approx. 16Gb) in the gvanno directory
   * `cd ~/gvanno`
   *  Download the [data bundle](https://drive.google.com/drive/folders/0B8aYD2TJ472mRUpFTEc4YzlTSUk) to `~/gvanno`
   * Decompress and untar the data bundle, e.g.: `tar -xvzf gvanno.bundle.v0.1.grch37.tgz`
3. Pull the gvanno Docker image from DockerHub:
   * `docker pull sigven/gvanno:latest`
4. Download the [gvanno pipeline script](https://github.com/sigven/gvanno/releases/download/v0.1/gvanno.sh) to `~/gvanno`

#### STEP 3: Run examples
