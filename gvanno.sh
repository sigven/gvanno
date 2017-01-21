#!/bin/bash

# -----------------------------------------
#                  Usage
# -----------------------------------------


#Set Script Name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`
CWD=$(pwd)

#Help function
function HELP {
  echo -e \\n"Help documentation for $SCRIPT"\\n
  echo -e "Usage: $SCRIPT GVANNO_DIRECTORY WORKING_DIRECTORY QUERY_VCF OUTPUT_VCF"\\n
  echo -e "Example: $SCRIPT ~/gvanno ~/gvanno sample1.vcf.gz sample1.annotated.vcf"\\n
  exit 1
}
SPLIT_MULTISAMPLE=0
#Check the number of arguments. If none are passed, print help and exit.
NUMARGS=$#
#echo $NUMARGS
if [ $NUMARGS -lt 4 ]; then
  HELP
fi

GVANNO_DIR=$1
WORKDIR=$2
CWD=$WORKDIR
QUERY_VCF=$3
QUERY_VCF_COMPLETE_PATH="$CWD/$3"
OUTPUT_VCF=$4
RAND_INDEX=$RANDOM
VEP_SUFFIX="/data/.vep"
VEPDB_DIR=$GVANNO_DIR$VEP_SUFFIX

#Make sure gvanno directory and VEP directory exists
if [ ! -d "$VEPDB_DIR" ] || [ ! -d "$GVANNO_DIR" ]; then
  echo -e \\n"ERROR: Directory $GVANNO_DIR or $VEPDB_DIR does not exist"\\n
  HELP
fi

if [ ! -f "$QUERY_VCF_COMPLETE_PATH" ] || [ ! -s "$QUERY_VCF_COMPLETE_PATH" ]; then
  echo -e \\n"ERROR: File $QUERY_VCF_COMPLETE_PATH	 does not exist or has zero size"\\n
  HELP
fi

BGZIP_SUFFIX=".gz"
if [[ "$QUERY_VCF" == *.vcf ]]
then
	VEP_OUTPUT_VCF=${QUERY_VCF/.vcf/.$RAND_INDEX.vep.vcf}
	VCFANNO_OUTPUT_VCF=${QUERY_VCF/.vcf/.$RAND_INDEX.vep.vcfanno.vcf}
	ANNOTATED_OUTPUT_VCF=${QUERY_VCF/.vcf/.$RAND_INDEX.vep.vcfanno.annotated.vcf}
fi
if [[ "$QUERY_VCF" == *.vcf.gz ]]
then
	VEP_OUTPUT_VCF=${QUERY_VCF/.vcf.gz/.$RAND_INDEX.vep.vcf}
	VCFANNO_OUTPUT_VCF=${QUERY_VCF/.vcf.gz/.$RAND_INDEX.vep.vcfanno.vcf}
	ANNOTATED_OUTPUT_VCF=${QUERY_VCF/.vcf.gz/.$RAND_INDEX.vep.vcfanno.annotated.vcf}
fi

VEP_OUTPUT_VCF_GZ=$VEP_OUTPUT_VCF$BGZIP_SUFFIX
VCFANNO_OUTPUT_VCF_GZ=$VCFANNO_OUTPUT_VCF$BGZIP_SUFFIX
ANNOTATED_OUTPUT_VCF_GZ=$ANNOTATED_OUTPUT_VCF$BGZIP_SUFFIX

VEP_OUTPUT_VCF_GZ_COMPLETE_PATH="$CWD/$VEP_OUTPUT_VCF_GZ"
VCFANNO_OUTPUT_VCF_GZ_COMPLETE_PATH="$CWD/$VCFANNO_OUTPUT_VCF_GZ"
ANNOTATED_OUTPUT_VCF_GZ_COMPLETE_PATH="$CWD/$ANNOTATED_OUTPUT_VCF_GZ"

VCFANNO_LOG=${VCFANNO_OUTPUT_VCF_GZ/.vcf.gz/.log}
OUTPUT_VCF_GZ=$OUTPUT_VCF$BGZIP_SUFFIX
OUTPUT_VCF_SAMPLE_POSTFIX=${OUTPUT_VCF_GZ/.vcf.gz/.vcf}
OUTPUT_TSV=${OUTPUT_VCF_GZ/.vcf.gz/.tsv}
VEP_OUTPUT_VCF_TMP=$(echo $VEP_OUTPUT_VCF | sed 's/.vcf/.tmp.vcf/')

CUST_USERID=""

current_date="`date +%Y-%m-%d`";
current_time="`date +%H:%M:%S`"
current_date_time="$current_date $current_time"
init_vep="$current_date_time - STEP 1: Basic variant annotation with Variant Effect Predictor (GENCODE, GRCh37)"
echo
echo $init_vep

#exit 1
docker run -i --rm $CUST_USERID -v=$VEPDB_DIR:/usr/local/share/vep/data -v=$CWD:/workdir -w=/workdir sigven/gvanno:latest sh -c "\
variant_effect_predictor.pl \
--input_file $QUERY_VCF \
--output_file $VEP_OUTPUT_VCF_TMP \
--vcf \
--check_ref \
--flag_pick_allele \
--force_overwrite \
--species homo_sapiens \
--assembly GRCh37 \
--offline \
--fork 4 \
--no_progress \
--variant_class \
--regulatory \
--domains \
--shift_hgvs 1 \
--hgvs \
--symbol \
--protein \
--ccds \
--uniprot \
--appris \
--biotype \
--canonical \
--gencode_basic \
--offline \
--cache \
--numbers \
--check_alleles \
--total_length \
--allele_number \
--no_escape \
--xref_refseq \
--dir /usr/local/share/vep/data \
--fasta /usr/local/share/vep/data/homo_sapiens/85_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"

docker run -i --rm $CUST_USERID -v=$VEPDB_DIR:/usr/local/share/vep/data -v=$CWD:/workdir -w=/workdir sigven/gvanno:latest sh -c "\
sed -r 's/\(p\.=\)|polyphen=2.2.2 |ClinVar=201507 |dbSNP=144 |ESP=20141103|sift=sift5.2.2 |COSMIC=71 |HGMD-PUBLIC=20152 //g' $VEP_OUTPUT_VCF_TMP > $VEP_OUTPUT_VCF"

docker run -i --rm $CUST_USERID -v=$VEPDB_DIR:/usr/local/share/vep/data -v=$CWD:/workdir -w=/workdir sigven/gvanno:latest sh -c "\
bgzip -f $VEP_OUTPUT_VCF"

docker run -i --rm $CUST_USERID -v=$VEPDB_DIR:/usr/local/share/vep/data -v=$CWD:/workdir -w=/workdir sigven/gvanno:latest sh -c "\
tabix -f -p vcf $VEP_OUTPUT_VCF_GZ"


rm -f $VEP_OUTPUT_VCF_TMP*

if [ ! -f "$VEP_OUTPUT_VCF_GZ_COMPLETE_PATH" ] || [ ! -s "$VEP_OUTPUT_VCF_GZ_COMPLETE_PATH" ]; then
   echo -e \\n"ERROR: File $VEP_OUTPUT_VCF_GZ_COMPLETE_PATH does not exist or has zero size"\\n
   HELP
fi

current_date="`date +%Y-%m-%d`";
current_time="`date +%H:%M:%S`"
current_date_time="$current_date $current_time"
init_vcfanno="$current_date_time - STEP 2:  Annotation for variant frequency and variant disease association using gvanno-vcfanno (ClinVar,dbSNP,dbNSFP,1000Genomes Project,ExAC,DoCM)"
echo
echo $init_vcfanno

docker run -i --rm $CUST_USERID -v=$GVANNO_DIR:/data -v=$CWD:/workdir -w=/workdir sigven/gvanno:latest sh -c "\
gvanno_vcfanno.py \
--dbsnp \
--dbnsfp \
--oneKG \
--docm \
--clinvar \
--exac \
$VEP_OUTPUT_VCF_GZ \
$VCFANNO_OUTPUT_VCF \
\$VCFANNO_DATA_DOCKER"

current_date="`date +%Y-%m-%d`";
current_time="`date +%H:%M:%S`"
current_date_time="$current_date $current_time"
echo "$current_date_time - Finished!"

if [ ! -f "$VCFANNO_OUTPUT_VCF_GZ_COMPLETE_PATH" ] || [ ! -s "$VCFANNO_OUTPUT_VCF_GZ_COMPLETE_PATH" ]; then
   echo -e \\n"ERROR: File $VCFANNO_OUTPUT_VCF_GZ_COMPLETE_PATH does not exist or has zero size"\\n
   HELP
fi

current_date="`date +%Y-%m-%d`";
current_time="`date +%H:%M:%S`"
current_date_time="$current_date $current_time"
init_summarise="$current_date_time - STEP 3: Gene annotations (phenotype associations, function) with gvanno-summarise"
echo
echo $init_summarise

docker run -i --rm $CUST_USERID -v=$GVANNO_DIR:/data -v=$CWD:/workdir -w=/workdir sigven/gvanno:latest sh -c "\
gvanno_summarise.py \
$VCFANNO_OUTPUT_VCF_GZ \
\$VCFANNO_DATA_DOCKER"

if [ ! -f "$ANNOTATED_OUTPUT_VCF_GZ_COMPLETE_PATH" ] || [ ! -s "$ANNOTATED_OUTPUT_VCF_GZ_COMPLETE_PATH" ]; then
   echo -e \\n"ERROR: Final VCF annotation file $ANNOTATED_OUTPUT_VCF_GZ_COMPLETE_PATH does not exist or has zero size"\\n
   HELP
fi

docker run -i --rm $CUST_USERID -v=$CWD:/workdir -w=/workdir sigven/gvanno:latest sh -c "\
bgzip -dc $ANNOTATED_OUTPUT_VCF_GZ > $OUTPUT_VCF"

docker run -i --rm $CUST_USERID -v=$CWD:/workdir -w=/workdir sigven/gvanno:latest sh -c "\
bgzip -f $OUTPUT_VCF"

docker run -i --rm $CUST_USERID -v=$CWD:/workdir -w=/workdir sigven/gvanno:latest sh -c "\
tabix -f -p vcf $OUTPUT_VCF_GZ"

#docker run -i --rm $CUST_USERID -v=$CWD:/workdir -w=/workdir sigven/gvanno:latest sh -c "\
#gvanno_vcf2tsv.py $OUTPUT_VCF_GZ $OUTPUT_TSV"

rm -f $VEP_OUTPUT_VCF_GZ_COMPLETE_PATH*
rm -f $VCFANNO_OUTPUT_VCF_GZ_COMPLETE_PATH*
rm -f $ANNOTATED_OUTPUT_VCF_GZ_COMPLETE_PATH*
rm -f "$CWD/$VCFANNO_LOG"
rm -f "$CWD/VEP_OUTPUT_VCF_TMP*"
