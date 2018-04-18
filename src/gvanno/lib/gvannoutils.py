#!/usr/bin/env python

import os,re,sys
import csv
import logging
import gzip
import toml

csv.field_size_limit(500 * 1024 * 1024)

def read_infotag_file(vcf_info_tags_tsv):
   """
   Function that reads a VCF info tag file that denotes annotation tags produced by PCGR.
   An example of the VCF info tag file is the following:
   
   tag	number	type	description
   Consequence	.	String	"Impact modifier for the consequence type (picked by VEP's --flag_pick_allele option)."
   
   A dictionary is returned, with the tag as the key, and the full dictionary record as the value
   """
   info_tag_xref = {} ##dictionary returned
   if not os.path.exists(vcf_info_tags_tsv):
      return info_tag_xref
   with open(vcf_info_tags_tsv, 'r') as tsvfile:
      reader = csv.DictReader(tsvfile, delimiter='\t')
      for rec in reader:
         if not rec['tag'] in info_tag_xref:
            info_tag_xref[rec['tag']] = rec
   
   return info_tag_xref


def gvanno_error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
   exit(1)

def read_config_options(configuration_file, gvanno_dir, genome_assembly, logger):
   
    ## read default options
   gvanno_config_options = {}
   gvanno_configuration_file_default = os.path.join(gvanno_dir,'data',str(genome_assembly),'gvanno_configuration_default.toml')
   if not os.path.exists(gvanno_configuration_file_default):
      err_msg = "Default gvanno configuration file " + str(gvanno_configuration_file_default) + " does not exist - exiting"
      gvanno_error_message(err_msg,logger)
   try:
      gvanno_config_options = toml.load(gvanno_configuration_file_default)
   except(IndexError,TypeError):
      err_msg = 'Configuration file ' + str(configuration_file) + ' is not formatted correctly'
      gvanno_error_message(err_msg, logger)

   ## override with options set by the users
   try:
      user_options = toml.load(configuration_file)
   except(IndexError,TypeError):
      err_msg = 'Configuration file ' + str(configuration_file) + ' is not formatted correctly'
      gvanno_error_message(err_msg, logger)
   
   
   boolean_tags = ['vep_skip_intergenic', 'vcf_validation']
   integer_tags = ['n_vcfanno_proc','n_vep_forks']
   for section in ['other']:
      if section in user_options:
         for t in boolean_tags:
            if t in user_options[section]:
               if not isinstance(user_options[section][t],bool):
                  err_msg = 'Configuration value ' + str(user_options[section][t]) + ' for ' + str(t) + ' cannot be parsed properly (expecting true/false)'
                  gvanno_error_message(err_msg, logger)
               gvanno_config_options[section][t] = int(user_options[section][t])
         for t in integer_tags:
            if t in user_options[section]:
               if not isinstance(user_options[section][t],int):
                  err_msg = 'Configuration value ' + str(user_options[section][t]) + ' for ' + str(t) + ' cannot be parsed properly (expecting integer)'
                  gvanno_error_message(err_msg, logger)
               gvanno_config_options[section][t] = user_options[section][t]
   
   return gvanno_config_options




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