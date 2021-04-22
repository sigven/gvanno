#!/usr/bin/env python

import csv
import re
import argparse
import os
import subprocess
import logging
import sys
import getpass
import platform
from argparse import RawTextHelpFormatter

GVANNO_VERSION = 'dev'
DB_VERSION = 'GVANNO_DB_VERSION = 20210422'
VEP_VERSION = '103'
GENCODE_VERSION = '37'
VEP_ASSEMBLY = "GRCh38"
DOCKER_IMAGE_VERSION = 'sigven/gvanno:' + str(GVANNO_VERSION)


def __main__():
   
   program_description = "gvanno - workflow for functional and clinical annotation of germline nucleotide variants"
   program_options = "   --query_vcf <QUERY_VCF>\n   --gvanno_dir <GVANNO_DIR>\n   --output_dir <OUTPUT_DIR>\n   --genome_assembly " + \
      "<grch37|grch38>\n   --sample_id <SAMPLE_ID>\n   --container <docker|singularity>"

   parser = argparse.ArgumentParser(description = program_description,
                                    formatter_class=RawTextHelpFormatter, usage="\n  %(prog)s -h [options]\n" + str(program_options))
   parser._action_groups.pop()
   required = parser.add_argument_group('Required arguments')
   optional_vep = parser.add_argument_group('VEP optional arguments')
   optional = parser.add_argument_group('Other optional arguments')

   optional.add_argument('--force_overwrite', action = "store_true", help='By default, the script will fail with an error if any ' + \
      'output file already exists.\nYou can force the overwrite of existing result files by using this flag, default: %(default)s')
   optional.add_argument('--version', action='version', version='%(prog)s ' + str(GVANNO_VERSION))
   optional.add_argument('--no_vcf_validate', action = "store_true",help="Skip validation of input VCF with Ensembl's vcf-validator, default: %(default)s")
   optional.add_argument('--docker-uid', dest = 'docker_user_id', help = 'Docker user ID. default is the host system user ID. If you are experiencing permission errors, try setting this up to root (`--docker-uid root`)')
   optional_vep.add_argument('--vep_regulatory', action='store_true', help = 'Enable Variant Effect Predictor (VEP) to look for overlap with regulatory regions (option --regulatory in VEP).')
   optional_vep.add_argument('--vep_lof_prediction', action = "store_true", help = "Predict loss-of-function variants with Loftee plugin " + \
      "in Variant Effect Predictor (VEP), default: %(default)s")
   optional_vep.add_argument('--vep_n_forks', default = 4, help="Number of forks for Variant Effect Predictor (VEP) processing, default: %(default)s")
   optional_vep.add_argument('--vep_buffer_size', default = 5000, help="Variant buffer size (variants read into memory simultaneously) " + \
      "for Variant Effect Predictor (VEP) processing\n- set lower to reduce memory usage, default: %(default)s")
   optional_vep.add_argument('--vep_pick_order', default = "canonical,appris,biotype,ccds,rank,tsl,length,mane", help="Comma-separated string " + \
      "of ordered transcript properties for primary variant pick in\nVariant Effect Predictor (VEP) processing, default: %(default)s")
   optional_vep.add_argument('--vep_skip_intergenic', action = "store_true", help="Skip intergenic variants in Variant Effect Predictor (VEP) processing, default: %(default)s")
   optional.add_argument('--vcfanno_n_processes', default = 4, help="Number of processes for vcfanno " + \
      "processing (see https://github.com/brentp/vcfanno#-p), default: %(default)s")
   required.add_argument('--query_vcf', help='VCF input file with germline query variants (SNVs/InDels).', required = True)
   required.add_argument('--gvanno_dir',help='Directory that contains the gvanno data bundle, e.g. ~/gvanno-' + str(GVANNO_VERSION), required = True)
   required.add_argument('--output_dir',help='Output directory', required = True)
   required.add_argument('--genome_assembly',choices = ['grch37','grch38'], help='Genome assembly build: grch37 or grch38', required = True)
   required.add_argument('--container', choices = ['docker', 'singularity'], action = "store",help="Run gvanno with docker or singularity")
   required.add_argument('--sample_id',help="Sample identifier - prefix for output files", required = True)
   
   args = parser.parse_args()
   arg_dict = vars(args)

   logger = getlogger('gvanno-check-workflow')

   ## Check that VEP pick criteria is formatted correctly
   if not arg_dict['vep_pick_order'] is None:
      values = str(arg_dict['vep_pick_order']).split(',')
      permitted_sources = ['canonical','appris','tsl','biotype','ccds','rank','length','mane']
      num_permitted_sources = 0
      for v in values:
         if v in permitted_sources:
            num_permitted_sources += 1
               
      if num_permitted_sources != 8:
         err_msg = "Option 'vep_pick_order' = " + str(arg_dict['vep_pick_order']) + " is formatted incorrectly, should be " + \
            "a comma-separated string of the following values: canonical,appris,tsl,biotype,ccds,rank,length,mane"
         gvanno_error_message(err_msg, logger)

   if arg_dict['container'] is None:
      err_msg = 'Please specify whether the gvanno workflow is running through Docker or Singularity (--container <docker|singularity>)'
      gvanno_error_message(err_msg, logger)

   logger = getlogger('gvanno-check-files')

   # check that script and Docker image version correspond
   if arg_dict['container'] == 'docker':
      check_docker_command = 'docker images -q ' + str(DOCKER_IMAGE_VERSION)
      output = subprocess.check_output(str(check_docker_command), stderr=subprocess.STDOUT, shell=True)
      
      if(len(output) == 0):
         err_msg = 'Docker image ' + str(DOCKER_IMAGE_VERSION) + ' does not exist, pull image from Dockerhub (docker pull ' + str(DOCKER_IMAGE_VERSION) + ')'
         gvanno_error_message(err_msg,logger)
   
   host_directories = verify_input_files(arg_dict, logger)

   run_gvanno(arg_dict, host_directories)


def gvanno_error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
   exit(0)

def verify_input_files(arg_dict, logger):
   """
   Function that checks the input files and directories provided by the user and checks for their existence
   """
 
   input_vcf_dir = "NA"
   input_conf_dir = "NA"
   db_dir = "NA"
   base_dir = "NA"
   output_dir_full = "NA"
   input_vcf_basename = "NA"
   input_conf_basename = "NA"
   
   ## check that query_vcf exist
   if arg_dict['query_vcf'] is None:
      err_msg = "Please specifiy a VCF input file (--query_vcf)"
      gvanno_error_message(err_msg,logger)
   
   ## check the existence of given output folder
   output_dir_full = os.path.abspath(arg_dict['output_dir'])
   if not os.path.isdir(output_dir_full):
      err_msg = "Output directory (" + str(output_dir_full) + ") does not exist"
      gvanno_error_message(err_msg,logger)
   
   ## check if input vcf exist
   if not arg_dict['query_vcf'] is None:
      if not os.path.exists(os.path.abspath(arg_dict['query_vcf'])):
         err_msg = "Input file (" + str(arg_dict['query_vcf']) + ") does not exist"
         gvanno_error_message(err_msg,logger)

      if not (os.path.abspath(arg_dict['query_vcf']).endswith('.vcf') or os.path.abspath(arg_dict['query_vcf']).endswith('.vcf.gz')):
         err_msg = "VCF input file (" + os.path.abspath(arg_dict['query_vcf']) + ") does not have the correct file extension (.vcf or .vcf.gz)"
         gvanno_error_message(err_msg,logger)

      ## check that tabix file exist if bgzipped files is given
      if os.path.abspath(arg_dict['query_vcf']).endswith('.vcf.gz'):
         tabix_file = arg_dict['query_vcf'] + '.tbi'
         if not os.path.exists(os.path.abspath(tabix_file)):
            err_msg = "Tabix file (i.e. '.gz.tbi') is not present for the bgzipped VCF input file (" + os.path.abspath(arg_dict['query_vcf']) + "). Please make sure your input VCF is properly compressed and indexed (bgzip + tabix)"
            gvanno_error_message(err_msg,logger)

      input_vcf_basename = os.path.basename(str(arg_dict['query_vcf']))
      input_vcf_dir = os.path.dirname(os.path.abspath(arg_dict['query_vcf']))

      ## if output vcf exist and overwrite not set
      output_vcf = os.path.join(str(output_dir_full),str(arg_dict['sample_id'])) + '_gvanno_' + str(arg_dict['genome_assembly']) + '.vcf.gz'
      if os.path.exists(output_vcf) and arg_dict['force_overwrite'] is False:
         err_msg = "Output files (e.g. " + str(output_vcf) + ") already exist - please specify different sample_id or add option --force_overwrite"
         gvanno_error_message(err_msg,logger)
   
   ## check the existence of base folder
   base_dir = os.path.abspath(arg_dict['gvanno_dir'])
   if not os.path.isdir(base_dir):
      err_msg = "Base directory (" + str(base_dir) + ") does not exist"
      gvanno_error_message(err_msg,logger)
   
   ## check the existence of data folder within the base folder
   db_dir = os.path.join(os.path.abspath(arg_dict['gvanno_dir']),'data')
   if not os.path.isdir(db_dir):
      err_msg = "Data directory (" + str(db_dir) + ") does not exist"
      gvanno_error_message(err_msg,logger)
   
   ## check the existence of specified assembly data folder within the base folder
   db_assembly_dir = os.path.join(os.path.abspath(arg_dict['gvanno_dir']),'data',arg_dict['genome_assembly'])
   if not os.path.isdir(db_assembly_dir):
      err_msg = "Data directory for the specified genome assembly (" + str(db_assembly_dir) + ") does not exist"
      gvanno_error_message(err_msg,logger)
   
   ## check the existence of RELEASE_NOTES (starting from 0.2.0)
   rel_notes_file = os.path.join(os.path.abspath(arg_dict['gvanno_dir']),'data',arg_dict['genome_assembly'],'RELEASE_NOTES')
   if not os.path.exists(rel_notes_file):
      err_msg = 'The gvanno data bundle is outdated - please download the latest data bundle (see github.com/sigven/gvanno for instructions)'
      gvanno_error_message(err_msg,logger)
   
   f_rel_not = open(rel_notes_file,'r')
   compliant_data_bundle = 0
   for line in f_rel_not:
      if DB_VERSION in line:
         compliant_data_bundle = 1
   
   f_rel_not.close()
         
   if compliant_data_bundle == 0:
      err_msg = 'The gvanno data bundle is not compliant with the software version - please download the latest software and data bundle (see https://github.com/sigven/gvanno for instructions)'
      gvanno_error_message(err_msg,logger)
   
   host_directories = {}
   host_directories['input_vcf_dir_host'] = input_vcf_dir
   host_directories['db_dir_host'] = db_assembly_dir
   host_directories['base_dir_host'] = base_dir
   host_directories['output_dir_host'] = output_dir_full
   host_directories['input_vcf_basename_host'] = input_vcf_basename
   
   return host_directories
   

def check_subprocess(command):
   try:
      output = subprocess.check_output(str(command), stderr=subprocess.STDOUT, shell=True)
      if len(output) > 0:
         print(str(output.decode()).rstrip())
   except subprocess.CalledProcessError as e:
      print(e.output.decode())
      exit(0)

def getlogger(logger_name):
   logger = logging.getLogger(logger_name)
   logger.setLevel(logging.DEBUG)

   # create console handler and set level to debug
   ch = logging.StreamHandler(sys.stdout)
   ch.setLevel(logging.DEBUG)

   # add ch to logger
   logger.addHandler(ch)
   
   # create formatter
   formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s", "20%y-%m-%d %H:%M:%S")
   
   #add formatter to ch
   ch.setFormatter(formatter)
   
   return logger

def run_gvanno(arg_dict, host_directories):
   """
   Main function to run the gvanno workflow using Docker
   """
   ## set basic Docker run commands
   output_vcf = 'None'
   output_pass_vcf = 'None'
   uid = ''
   docker_user_id = arg_dict['docker_user_id']

   global GENCODE_VERSION, VEP_ASSEMBLY
   if arg_dict['genome_assembly'] == 'grch37':
      GENCODE_VERSION = 'release 19'
      VEP_ASSEMBLY = 'GRCh37'

   logger = getlogger('gvanno-get-OS')
   if docker_user_id:
      uid = docker_user_id
   elif platform.system() == 'Linux' or platform.system() == 'Darwin' or sys.platform == 'darwin' or sys.platform == 'linux2' or sys.platform == 'linux':
      uid = os.getuid()
   else:
      if platform.system() == 'Windows' or sys.platform == 'win32' or sys.platform == 'cygwin':
         uid = getpass.getuser()
   
   if uid == '':
      logger.warning('Was not able to get user id/username for logged-in user on the underlying platform (platform.system(): ' + str(platform.system()) + ', sys.platform: ' + str(sys.platform) + '), now running gvanno as root')
      uid = 'root'
      
   vepdb_dir_host = os.path.join(str(host_directories['db_dir_host']),'.vep')
   vcf_validation = 1
   if arg_dict['no_vcf_validate']:
      vcf_validation = 0
   data_dir = '/data'
   output_dir = '/workdir/output'
   vep_dir = '/usr/local/share/vep/data'
   input_vcf_docker = 'None'
   
   if host_directories['input_vcf_basename_host'] != 'NA':
      input_vcf_docker = '/workdir/input_vcf/' + str(host_directories['input_vcf_basename_host'])
   
   vep_volume_mapping = str(vepdb_dir_host) + ":/usr/local/share/vep/data"
   databundle_volume_mapping = str(host_directories['base_dir_host']) + ":/data"
   input_vcf_volume_mapping = str(host_directories['input_vcf_dir_host']) + ":/workdir/input_vcf"
   output_volume_mapping = str(host_directories['output_dir_host']) + ":/workdir/output"

   if arg_dict['container'] == 'docker':
      container_command_run1 = "docker run --rm -t -u " + str(uid) + " -v=" +  str(databundle_volume_mapping) + " -v=" + str(vep_volume_mapping) + " -v=" + str(output_volume_mapping)
   elif arg_dict['container'] == 'singularity':
      container_command_run1 = "singularity exec " + " -B " +  str(databundle_volume_mapping) + " -B " + str(vep_volume_mapping) + " -B " + str(output_volume_mapping)

   if host_directories['input_vcf_dir_host'] != 'NA' and arg_dict['container'] == 'docker':
      container_command_run1 = container_command_run1  + " -v=" + str(input_vcf_volume_mapping)
   elif host_directories['input_vcf_dir_host'] != 'NA' and arg_dict['container'] == 'singularity':
      container_command_run1 = container_command_run1  + " -B " + str(input_vcf_volume_mapping)

   if arg_dict['container'] == 'docker':
      container_command_run1 = container_command_run1 + " -w=/workdir/output " + str(DOCKER_IMAGE_VERSION) + " sh -c \""
   elif arg_dict['container'] == 'singularity':
      container_command_run1 = container_command_run1 + " -W /workdir/output " + 'src/gvanno.simg' + " sh -c \""

   if arg_dict['container'] == 'docker':
      container_command_run2 = "docker run --rm -t -u " + str(uid) + " -v=" +  str(databundle_volume_mapping) + " -v=" + str(output_volume_mapping)
      container_command_run2 = container_command_run2 + " -w=/workdir/output " + str(DOCKER_IMAGE_VERSION) + " sh -c \""
      docker_command_run_end = '\"'
   elif arg_dict['container'] == 'singularity':
      container_command_run2 = "singularity exec " + " -B " +  str(databundle_volume_mapping) + " -B " + str(output_volume_mapping)
      container_command_run2 = container_command_run2 + " -W /workdir/output " + 'src/gvanno.simg' + " sh -c \""
      docker_command_run_end = '\"'

   ## GVANNO|start - Log key information about sample, options and assembly
   logger = getlogger("gvanno-start")
   logger.info("--- Germline variant annotation (gvanno) workflow ----")
   logger.info("Sample name: " + str(arg_dict['sample_id']))
   logger.info("Genome assembly: " + str(arg_dict['genome_assembly']))
   print()

   ## GVANNO|validate - verify input file (contents/format)
   logger = getlogger('gvanno-validate-input')
   logger.info("STEP 0: Validate input data")
   vcf_validate_command = str(container_command_run1) + "gvanno_validate_input.py " + str(data_dir) + " " + str(input_vcf_docker) + " " + \
      str(vcf_validation) + " "  + str(arg_dict['genome_assembly']) + docker_command_run_end

   check_subprocess(vcf_validate_command)
   logger.info('Finished')
   
   if not input_vcf_docker == 'None':
      
      ## Define input, output and temporary file names
      output_vcf = os.path.join(output_dir, str(arg_dict['sample_id']) + '_gvanno_' + str(arg_dict['genome_assembly']) + '.vcf.gz')
      output_tsv = os.path.join(output_dir, str(arg_dict['sample_id']) + '_gvanno_'  + str(arg_dict['genome_assembly']) + '.tsv')
      output_pass_vcf = os.path.join(output_dir, str(arg_dict['sample_id']) + '_gvanno_pass_' + str(arg_dict['genome_assembly']) + '.vcf.gz')
      output_pass_tsv = os.path.join(output_dir, str(arg_dict['sample_id']) + '_gvanno_pass_' + str(arg_dict['genome_assembly']) + '.tsv')
      input_vcf_gvanno_ready = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.gvanno_ready.vcf.gz',host_directories['input_vcf_basename_host']))
      vep_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.vep.vcf',input_vcf_gvanno_ready)
      vep_vcfanno_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.vep.vcfanno.vcf',input_vcf_gvanno_ready)
      vep_vcfanno_annotated_vcf = re.sub(r'\.vcfanno','.vcfanno.annotated',vep_vcfanno_vcf) + '.gz'
      vep_vcfanno_annotated_pass_vcf = re.sub(r'\.vcfanno','.vcfanno.annotated.pass',vep_vcfanno_vcf) + '.gz'
      
      ## Path for human genome assembly and human ancestor (FASTA)
      fasta_assembly = os.path.join(vep_dir, "homo_sapiens", str(VEP_VERSION) + "_" + str(VEP_ASSEMBLY), "Homo_sapiens." + str(VEP_ASSEMBLY) + ".dna.primary_assembly.fa.gz")
      ancestor_assembly = os.path.join(vep_dir, "homo_sapiens", str(VEP_VERSION) + "_" + str(VEP_ASSEMBLY), "human_ancestor.fa.gz")

      ## List all VEP flags used when calling VEP
      loftee_dir = '/opt/vep/src/ensembl-vep/modules'
      plugins_in_use = "NearestExonJB"
      vep_flags = "--hgvs --dont_skip --failed 1 --af --af_1kg --af_gnomad --variant_class --domains --symbol --protein --ccds " + \
         "--uniprot --appris --biotype --canonical --gencode_basic --mane --cache --numbers --total_length --allele_number --no_escape " + \
         "--xref_refseq --plugin NearestExonJB,max_range=50000"
      vep_options = "--vcf --quiet --check_ref --flag_pick_allele --pick_order " + str(arg_dict['vep_pick_order']) + \
         " --force_overwrite --species homo_sapiens --assembly " + str(VEP_ASSEMBLY) + " --offline --fork " + \
         str(arg_dict['vep_n_forks']) + " " + str(vep_flags) + " --dir /usr/local/share/vep/data"
      if arg_dict['vep_skip_intergenic'] == 1:
         vep_options = vep_options + " --no_intergenic"
      if arg_dict['vep_regulatory'] == 1:
         vep_options = vep_options + " --regulatory"
      if arg_dict['vep_lof_prediction'] == 1:
         plugins_in_use = plugins_in_use + ", LoF"
         vep_options += " --plugin LoF,loftee_path:" + loftee_dir + ",human_ancestor_fa:" + str(ancestor_assembly)  + ",use_gerp_end_trunc:0 --dir_plugins " + loftee_dir

      ## Compose full VEP command
      vep_main_command = str(container_command_run1) + "vep --input_file " + str(input_vcf_gvanno_ready) + " --output_file " + str(vep_vcf) + \
         " " + str(vep_options) + " --buffer_size " + str(arg_dict['vep_buffer_size']) + " --fasta " + str(fasta_assembly) + docker_command_run_end
      vep_bgzip_command = container_command_run1 + "bgzip -f -c " + str(vep_vcf) + " > " + str(vep_vcf) + ".gz" + docker_command_run_end
      vep_tabix_command = str(container_command_run1) + "tabix -f -p vcf " + str(vep_vcf) + ".gz" + docker_command_run_end

      ## GVANNO|VEP - run consequence annotation with Variant Effect Predictor
      logger = getlogger('gvanno-vep')   
      print()
      logger.info("STEP 1: Basic variant annotation with Variant Effect Predictor (" + str(VEP_VERSION) + ", GENCODE " + str(GENCODE_VERSION) + ", " + str(arg_dict['genome_assembly']) + ")")
      logger.info("VEP configuration - one primary consequence block pr. alternative allele (--flack_pick_allele)")
      logger.info("VEP configuration - transcript pick order: " + str(arg_dict['vep_pick_order']))
      logger.info("VEP configuration - transcript pick order: See more at https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options")
      logger.info("VEP configuration - buffer size: " + str(arg_dict['vep_buffer_size']))
      logger.info("VEP configuration - skip intergenic: " + str(arg_dict['vep_skip_intergenic']))
      logger.info("VEP configuration - look for overlap with regulatory regions: " + str(arg_dict['vep_regulatory']))
      logger.info("VEP configuration - number of forks: " + str(arg_dict['vep_n_forks']))
      logger.info("VEP configuration - loss-of-function prediction: " + str(arg_dict['vep_lof_prediction']))
      logger.info("VEP configuration - plugins in use: " + str(plugins_in_use))

      check_subprocess(vep_main_command)
      check_subprocess(vep_bgzip_command)
      check_subprocess(vep_tabix_command)
      logger.info("Finished")

      ## GVANNO|vcfanno - annotate VCF against a number of variant annotation resources
      print()
      logger = getlogger('gvanno-vcfanno')
      logger.info("STEP 2: Clinical/functional variant annotations with gvanno-vcfanno (ClinVar, ncER, dbNSFP, GWAS catalog, UniProtKB, cancerhotspots.org)")
      logger.info('vcfanno configuration - number of processes (-p): ' + str(arg_dict['vcfanno_n_processes']))
      gvanno_vcfanno_command = str(container_command_run2) + "gvanno_vcfanno.py --num_processes "  + str(arg_dict['vcfanno_n_processes']) + \
         " --dbnsfp --clinvar --ncer --uniprot --gvanno_xref --gwas --cancer_hotspots " + str(vep_vcf) + ".gz " + str(vep_vcfanno_vcf) + \
         " " + os.path.join(data_dir, "data", str(arg_dict['genome_assembly'])) + docker_command_run_end
      check_subprocess(gvanno_vcfanno_command)
      logger.info("Finished")

      ## GVANNO|summarise - expand annotations in VEP and vcfanno-annotated VCF file
      print()
      logger = getlogger("gvanno-summarise")
      logger.info("STEP 3: Summarise gene and variant annotations with gvanno-summarise")
      gvanno_summarise_command = str(container_command_run2) + "gvanno_summarise.py " + str(vep_vcfanno_vcf) + ".gz " + \
         os.path.join(data_dir, "data", str(arg_dict['genome_assembly'])) + " " + str(int(arg_dict['vep_lof_prediction'])) + docker_command_run_end
      check_subprocess(gvanno_summarise_command)
      logger.info("Finished")
      
      ## GVANNO|clean - move output files and clean up temporary files
      create_output_vcf_command1 = str(container_command_run2) + 'mv ' + str(vep_vcfanno_annotated_vcf) + ' ' + str(output_vcf) + "\""
      create_output_vcf_command2 = str(container_command_run2) + 'mv ' + str(vep_vcfanno_annotated_vcf) + '.tbi ' + str(output_vcf) + '.tbi' + "\""
      create_output_vcf_command3 = str(container_command_run2) + 'mv ' + str(vep_vcfanno_annotated_pass_vcf) + ' ' + str(output_pass_vcf) + "\""
      create_output_vcf_command4 = str(container_command_run2) + 'mv ' + str(vep_vcfanno_annotated_pass_vcf) + '.tbi ' + str(output_pass_vcf) + '.tbi' + "\""
      clean_command = str(container_command_run2) + 'rm -f ' + str(vep_vcf) + '* ' + str(vep_vcfanno_annotated_vcf) + ' ' + \
          str(vep_vcfanno_annotated_pass_vcf) + '* ' + str(vep_vcfanno_vcf) + '* ' +  str(input_vcf_gvanno_ready) + "* "  + docker_command_run_end
      check_subprocess(create_output_vcf_command1)
      check_subprocess(create_output_vcf_command2)
      check_subprocess(create_output_vcf_command3)
      check_subprocess(create_output_vcf_command4)
      check_subprocess(clean_command)
      
      print()
      ## GVANNO|vcf2tsv - convert VCF to TSV with https://github.com/sigven/vcf2tsv
      logger = getlogger("gvanno-vcf2tsv")
      logger.info("STEP 4: Converting VCF to TSV with https://github.com/sigven/vcf2tsv")
      gvanno_vcf2tsv_command_pass = str(container_command_run2) + "vcf2tsv.py " + str(output_pass_vcf) + " --compress " + str(output_pass_tsv) + docker_command_run_end
      gvanno_vcf2tsv_command_all = str(container_command_run2) + "vcf2tsv.py " + str(output_vcf) + " --compress --keep_rejected " + str(output_tsv) + docker_command_run_end
      logger.info("Conversion of VCF variant data to records of tab-separated values - PASS variants only")
      check_subprocess(gvanno_vcf2tsv_command_pass)
      logger.info("Conversion of VCF variant data to records of tab-separated values - PASS and non-PASS variants")
      check_subprocess(gvanno_vcf2tsv_command_all)
      logger.info("Finished")
      
      #return
  
   print
   
   
if __name__=="__main__": __main__()
