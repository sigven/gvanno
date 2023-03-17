#!/usr/bin/env python

import argparse
import re
import os
import subprocess
import logging
import sys
import locale
#import wget
import urllib.request as urllib2
from argparse import RawTextHelpFormatter

GVANNO_VERSION = '1.6.0'
REFDATA_VERSION = '20230310'
ENSEMBL_VERSION = '109'
GENCODE_VERSION = 'v43'
VEP_ASSEMBLY = "GRCh38"
HOST_GVANNO_REFDATA_URL = "http://insilico.hpc.uio.no/pcgr/gvanno/"

def __main__():
   
   program_description = "download_gvanno_refdata - download reference datasets for the gvanno variant annotation workflow"
   program_options = "   --download_dir <DOWNLOAD_DIR>\n   --genome_assembly " + \
      "<grch37|grch38>\n"

   parser = argparse.ArgumentParser(description = program_description,
                                    formatter_class=RawTextHelpFormatter, usage="\n  %(prog)s -h [options]\n" + str(program_options))
   parser._action_groups.pop()
   required = parser.add_argument_group('Required arguments')
   optional = parser.add_argument_group('Other optional arguments')

   optional.add_argument('--force_overwrite', action = "store_true", help='By default, the script will fail with an error if any ' + \
      'download directory already exist.\nYou can force the overwrite of existing download directory by using this flag, default: %(default)s')
   optional.add_argument('--version', action='version', version='%(prog)s ' + str(GVANNO_VERSION))
   optional.add_argument("--debug", action="store_true", help="Print full commands to log and do not delete intermediate files with warnings etc.")
   required.add_argument('--download_dir',help='Destination directory for downloaded reference data', required = True)
   required.add_argument('--genome_assembly',choices = ['grch37','grch38'], help='Choose build-specific reference data for download: grch37 or grch38', required = True)

   
   args = parser.parse_args()
   arg_dict = vars(args)

   logger = getlogger('check-download-dir')

   download_dir_full = os.path.abspath(arg_dict['download_dir'])
   if not os.path.isdir(download_dir_full):
      err_msg = "Download destination directory (" + str(download_dir_full) + ") does not exist"
      gvanno_error_message(err_msg,logger)
   
   ## check that 'data' does not exist in download_dir_full
   db_dir = os.path.join(os.path.abspath(arg_dict['download_dir']),'data')
   arg_dict['db_assembly_dir'] = os.path.join(os.path.abspath(arg_dict['download_dir']),'data',arg_dict['genome_assembly'])
   arg_dict['vep_assembly_dir'] = os.path.join(os.path.abspath(arg_dict['download_dir']),'data',arg_dict['genome_assembly'], '.vep')

   if os.path.isdir(db_dir) and arg_dict['force_overwrite'] is False:
      err_msg = "Data directory (" + str(db_dir) + ") exists (force_overwrite = False)"
      gvanno_error_message(err_msg,logger)
   
   if not os.path.isdir(db_dir):
      os.mkdir(db_dir)
   if not os.path.isdir(arg_dict['db_assembly_dir']):
      os.mkdir(arg_dict['db_assembly_dir'])
      os.mkdir(arg_dict['vep_assembly_dir'])


   download_gvanno_ref_data(arg_dict = arg_dict)


def gvanno_error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
   exit(0)

def check_subprocess(command, logger, debug = False):
   if debug:
      logger.info(command)
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

def get_url_num_bytes(url,  logger, verbose=True):
   """
   Get number of bytes of file at a URL. None if not reported.
   """
   # urllib2 could be a bit more intelligent in guessing what I mean:
   if not re.match('^[a-zA-Z]*:', url):
      if os.path.exists(url):
         url = 'file:' + url
      else:
         url = 'http://' + url
   try:
      class HeadRequest(urllib2.Request):
         def get_method(self):
            return "HEAD"
      response = urllib2.urlopen(HeadRequest(
               url, None, {'User-Agent': 'Mozilla/5.0'}))
       
      num_bytes = response.info().get('content-length')
      if verbose:
         # urllib does redirections for me. Report if it's been clever:
         if response.geturl() != url:
            logger.info('Reporting size of file at: ' + response.geturl())
   except:
      logger.info('Failed to connect to URL ' + str(url))
      num_bytes = None
   if num_bytes is not None:
      num_bytes = int(num_bytes)
   return num_bytes

def pretty_print(num_bytes, logger):
   """
   Output number of bytes according to locale and with IEC binary prefixes
   """
   if num_bytes is None:
      logger.info('File size unavailable.')
      return
   KiB = 1024
   MiB = KiB * KiB
   GiB = KiB * MiB
   TiB = KiB * GiB
   PiB = KiB * TiB
   EiB = KiB * PiB
   ZiB = KiB * EiB
   YiB = KiB * ZiB
   locale.setlocale(locale.LC_ALL, '')
   output = locale.format_string("%d", num_bytes, grouping=True) + ' bytes'
   if num_bytes > YiB:
      output += ' (%.3g YiB)' % (num_bytes / YiB)
   elif num_bytes > ZiB:
      output += ' (%.3g ZiB)' % (num_bytes / ZiB)
   elif num_bytes > EiB:
      output += ' (%.3g EiB)' % (num_bytes / EiB)
   elif num_bytes > PiB:
      output += ' (%.3g PiB)' % (num_bytes / PiB)
   elif num_bytes > TiB:
      output += ' (%.3g TiB)' % (num_bytes / TiB)
   elif num_bytes > GiB:
      output += ' (%.3g GiB)' % (num_bytes / GiB)
   elif num_bytes > MiB:
      output += ' (%.3g MiB)' % (num_bytes / MiB)
   elif num_bytes > KiB:
      output += ' (%.3g KiB)' % (num_bytes / KiB)
   
   return(output)
   #logger.info(output)

def download_gvanno_ref_data(arg_dict):
   """
   Main function to run the gvanno workflow using Docker
   """
   global GENCODE_VERSION, VEP_ASSEMBLY
   if arg_dict['genome_assembly'] == 'grch37':
      GENCODE_VERSION = 'v19'
      VEP_ASSEMBLY = 'GRCh37'

   vep_assembly_dir = os.path.join(os.path.abspath(arg_dict['download_dir']),'data',arg_dict['genome_assembly'], '.vep')
   
   datasets = {}
   for db in ['vep_cache','vep_fasta','gvanno_custom']:
      datasets[db] = {}
      datasets[db]['remote_url'] = 'NA'
      datasets[db]['local_path'] = 'NA'

   ## Remote URL - pre-indexed VEP cache
   datasets['vep_cache']['remote_url'] = (
      f"ftp://ftp.ensembl.org/pub/current_variation/indexed_vep_cache/"
      f"homo_sapiens_vep_{ENSEMBL_VERSION}_{VEP_ASSEMBLY}.tar.gz"
      )
   
   ## Remote URL - Human genome FASTA
   datasets['vep_fasta']['remote_url'] = (
      f"http://ftp.ensembl.org/pub/release-{ENSEMBL_VERSION}/fasta/homo_sapiens/dna/Homo_sapiens."
      f"{VEP_ASSEMBLY}.dna.primary_assembly.fa.gz"
      )
   
   logger = getlogger('download-vep-cache')
   datasets['vep_cache']['local_path'] = os.path.join(arg_dict['vep_assembly_dir'], f"homo_sapiens_vep_{ENSEMBL_VERSION}_{VEP_ASSEMBLY}.tar.gz")
   datasets['vep_fasta']['local_path'] = os.path.join(arg_dict['vep_assembly_dir'], "homo_sapiens", f"{ENSEMBL_VERSION}_{VEP_ASSEMBLY}", f"Homo_sapiens.{VEP_ASSEMBLY}.dna.primary_assembly.fa.gz")
   datasets['vep_fasta']['local_path_uncompressed'] = re.sub(r'.gz','',datasets['vep_fasta']['local_path'])

   vep_cache_bytes_remote = get_url_num_bytes(url = datasets['vep_cache']['remote_url'], logger = logger)

   if VEP_ASSEMBLY == 'GRCh37':
      datasets['vep_fasta']['remote_url'] = (
      f"http://ftp.ensembl.org/pub/grch37/release-{ENSEMBL_VERSION}/fasta/homo_sapiens/dna/Homo_sapiens."
      f"{VEP_ASSEMBLY}.dna.primary_assembly.fa.gz"
      )
   
  
   logger.info('VEP cache - remote target file ' + str(datasets['vep_cache']['remote_url']))
   logger.info('VEP cache - size: ' + pretty_print(vep_cache_bytes_remote, logger = logger))
   logger.info('VEP cache - local destination file: ' + str(datasets['vep_cache']['local_path']))

   if os.path.exists(datasets['vep_cache']['local_path']):
      if os.path.getsize(datasets['vep_cache']['local_path']) == vep_cache_bytes_remote:
         logger.info('VEP cache already downloaded')
      else:
         logger.info('VEP cache - download in progress - this can take a while ...  ')
         #urllib2.urlretrieve(datasets['vep_cache']['remote_url'], datasets['vep_cache']['local_path'])
   else:
      logger.info('VEP cache - download in progress - this can take a while ...  ')
      #urllib2.urlretrieve(datasets['vep_cache']['remote_url'], datasets['vep_cache']['local_path'])
   
   logger.info('VEP cache - unzip and untar')

   wdir = os.getcwd()
   
   os.chdir(vep_assembly_dir)
   
   command_unzip_cache = f"gzip -dc {datasets['vep_cache']['local_path']} | tar xvf -"
   check_subprocess(command = command_unzip_cache, logger = logger)

   ## change back to working directory
   os.chdir(wdir)

   logger = getlogger('download-vep-fasta')
   fasta_cache_bytes = get_url_num_bytes(url = datasets['vep_fasta']['remote_url'], logger = logger)
   logger.info('VEP FASTA - remote target file: ' + str(datasets['vep_fasta']['remote_url']))
   logger.info('VEP FASTA - size: ' + pretty_print(fasta_cache_bytes, logger = logger))
   logger.info('VEP FASTA - local destination file: ' + str(datasets['vep_fasta']['local_path']))
   logger.info('VEP FASTA - download in progress - this can take a while ...  ')

   #urllib2.urlretrieve(datasets['vep_fasta']['remote_url'], datasets['vep_fasta']['local_path'])
   logger.info('VEP FASTA - unzip + bgzip')
   command_unzip_fasta = f"gzip -d {datasets['vep_fasta']['local_path']}"
   check_subprocess(command = command_unzip_fasta, logger = logger)
   command_bgzip_fasta = f"bgzip {datasets['vep_fasta']['local_path_uncompressed']}"
   check_subprocess(command = command_bgzip_fasta, logger = logger)

   datasets['gvanno_custom']['remote_url'] = f'{HOST_GVANNO_REFDATA_URL}gvanno.databundle.{arg_dict["genome_assembly"]}.{REFDATA_VERSION}.tgz'
   datasets['gvanno_custom']['local_path'] = os.path.join(arg_dict['download_dir'], f'gvanno.databundle.{arg_dict["genome_assembly"]}.{REFDATA_VERSION}.tgz')

   logger = getlogger('download-gvanno-custom')
   custom_cache_bytes = get_url_num_bytes(url = datasets['gvanno_custom']['remote_url'], logger = logger)
   logger.info("Downloading custom gvanno variant datasets:  Clinvar / dbNSFP / ncER / cancerhotspots ++")
   logger.info('Custom gvanno datasets - remote target file ' + str(datasets['gvanno_custom']['remote_url']))
   logger.info('Custom gvanno datasets - size: ' + pretty_print(custom_cache_bytes, logger = logger))

   urllib2.urlretrieve(datasets['gvanno_custom']['remote_url'], datasets['custom_gvanno']['local_path'])
   logger.info('Custom gvanno datasets - unzip and untar')
   
   wdir = os.getcwd()
   os.chdir(arg_dict['download_dir'])
   
   command_unzip_cache = f"gzip -dc {datasets['gvanno_custom']['local_path']} | tar xvf -"

   logger.info('Finished')
   
   
if __name__=="__main__": __main__()
