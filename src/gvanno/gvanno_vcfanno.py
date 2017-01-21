#!/usr/bin/env python

import argparse
import cyvcf
import vcfutils
import random
import gvanno
import os
import re
import sys

logger = gvanno.getlogger('gvanno-vcfanno')


def __main__():
   parser = argparse.ArgumentParser(description='Run brentp/vcfanno - annotate a VCF file against multiple VCF files in parallel', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('query_vcf', help='Bgzipped input VCF file with query variants (SNVs/InDels)')
   parser.add_argument('out_vcf', help='Output VCF file with appended annotations from multiple VCF files')
   parser.add_argument('gvanno_dir', help='gvanno base directory')
   parser.add_argument("--exac",action = "store_true", help="Annotate VCF with annotations from Exome Aggregation Consortium (release 0.3.1)")
   parser.add_argument("--docm",action = "store_true", help="Annotate VCF with annotations from Database of Curated Mutations (DoCM 3.2, May 2016)")
   parser.add_argument("--clinvar",action = "store_true", help="Annotate VCF with annotations from ClinVar (October 2016)")
   parser.add_argument("--dbsnp",action = "store_true", help="Annotate VCF with annotations from database of short genetic variations (dbSNP b147)")
   parser.add_argument("--dbnsfp",action = "store_true", help="Annotate VCF with annotations from database of non-synonymous functional predictions (dbNSFP v3.2)")
   parser.add_argument("--oneKG",action = "store_true", help="Annotate VCF with annotations from the 1000 Genome Project (1000GProject wgs phase 3)")

      
   args = parser.parse_args()

   query_info_tags = get_vcf_info_tags(args.query_vcf)
   vcfheader_file = args.out_vcf + '.tmp.' + str(random.randrange(0,10000000)) + '.header.txt'
   conf_fname = args.out_vcf + '.tmp.conf.toml'
   print_vcf_header(args.query_vcf, vcfheader_file, chromline_only = False)
   run_vcfanno(args.query_vcf, query_info_tags, vcfheader_file, args.gvanno_dir, conf_fname, args.out_vcf, args.exac, args.docm, args.clinvar, args.dbsnp, args.dbnsfp, args.oneKG)

def run_vcfanno(query_vcf, query_info_tags, vcfheader_file, gvanno_directory, conf_fname, output_vcf, exac, docm, clinvar, dbsnp, dbnsfp, oneKG):
   
   gvanno_db_directory = gvanno_directory + '/data'
   if exac is True:
      exac_tags = ["AMR_AF_EXAC","AFR_AF_EXAC","NFE_AF_EXAC","FIN_AF_EXAC","OTH_AF_EXAC","GLOBAL_AF_EXAC","EAS_AF_EXAC","SAS_AF_EXAC"]
      for t in exac_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the ExAC VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("exac",gvanno_db_directory, conf_fname)
      append_to_vcf_header(gvanno_db_directory, "exac", vcfheader_file)
   if docm is True:
      docm_tags = ["DOCM_DISEASE","DOCM_PMID"]
      for t in docm_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the DoCM VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("docm",gvanno_db_directory, conf_fname)
      append_to_vcf_header(gvanno_db_directory, "docm", vcfheader_file)

   if clinvar is True:
      clinvar_tags = ["CLINVAR_MSID","CLINVAR_PMIDS","CLINVAR_SIG","CLINVAR_VARIANT_ORIGIN"]
      for t in clinvar_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the ClinVar VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("clinvar",gvanno_db_directory, conf_fname)
      append_to_vcf_header(gvanno_db_directory, "clinvar", vcfheader_file)
   if dbsnp is True:
      dbsnp_tags = ["GWAS_CATALOG_PMID","GWAS_CATALOG_TRAIT_URI","DBSNPRSID", "DBSNPBUILDID", "DBSNP_VALIDATION","DBSNP_MAPPINGSTATUS","DBSNP_SUBMISSIONS"]
      for t in dbsnp_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the dbSNP VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("dbsnp",gvanno_db_directory, conf_fname)
      append_to_vcf_header(gvanno_db_directory, "dbsnp", vcfheader_file)
   if dbnsfp is True:
      dbnsfp_tags = ["DBNSFP"]
      for t in dbnsfp_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the dbNSFP VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("dbnsfp",gvanno_db_directory, conf_fname)
      append_to_vcf_header(gvanno_db_directory, "dbnsfp", vcfheader_file)
   if oneKG is True:
      oneKG_tags = ["EAS_AF_1KG","EUR_AF_1KG","AMR_AF_1KG","AFR_AF_1KG","SAS_AF_1KG","GLOBAL_AF_1KG"]
      for t in oneKG_tags:
         if query_info_tags.has_key(t):
            logger.warning("Query VCF has INFO tag " + str(t) + ' - this is also present in the 1000Genomes Project VCF annotation file. This tag will be overwritten if not renamed in the query VCF')
      append_to_conf_file("oneKG",gvanno_db_directory, conf_fname)
      append_to_vcf_header(gvanno_db_directory, "oneKG", vcfheader_file)

   out_vcf_vcfanno_unsorted1 = output_vcf + '.tmp.unsorted.1'
   out_vcf_vcfanno_unsorted2 = output_vcf + '.tmp.unsorted.2'
   out_vcf_vcfanno_sorted = output_vcf + '.tmp.sorted.1'
   query_prefix = re.sub('\.vcf.gz$','',query_vcf)
   print_vcf_header(query_vcf, vcfheader_file, chromline_only = True)
   command1 = "vcfanno -p=4 " + str(conf_fname) + " " + str(query_vcf) + " > " + str(out_vcf_vcfanno_unsorted1) + " 2> " + str(query_prefix) + '.vcfanno.log'
   os.system(command1)
   
   os.system('cat ' + str(vcfheader_file) + ' > ' + str(output_vcf))
   os.system('cat ' + str(out_vcf_vcfanno_unsorted1) + ' | grep -v \'^#\' >> ' + str(output_vcf))
   os.system('rm -f ' + str(output_vcf) + '.tmp*')
   os.system('bgzip ' + str(output_vcf))
   os.system('tabix -p vcf ' + str(output_vcf) + '.gz')
   return 0
   
def append_to_vcf_header(gvanno_db_directory, dbsource, vcfheader_file):
   
   vcf_info_tags_file = str(gvanno_db_directory) + '/' + str(dbsource) + '/' + str(dbsource) + '.vcfanno.vcf_info_tags.txt'
   os.system('cat ' + str(vcf_info_tags_file) + ' >> ' + str(vcfheader_file))


def append_to_conf_file(dbsource, gvanno_db_directory, conf_fname):
   fh = open(conf_fname,'a')
   fh.write('[[annotation]]\n')
   fh.write('file="' + str(gvanno_db_directory) + '/' + str(dbsource) + '/' + str(dbsource) + '.vcf.gz"\n')
   if dbsource == 'dbsnp':
      fh.write('fields = ["GWAS_CATALOG_PMID","GWAS_CATALOG_TRAIT_URI","DBSNPRSID", "DBSNPBUILDID", "DBSNP_VALIDATION","DBSNP_MAPPINGSTATUS","DBSNP_SUBMISSIONS"]\n')
      fh.write('ops=["concat","concat", "concat", "concat", "concat", "concat","concat"]\n\n')
   if dbsource == 'dbnsfp':
      
      fh.write('fields = ["DBNSFP"]\n')
      fh.write('ops=["concat"]\n\n')
      
   if dbsource == 'exac':
      fh.write('fields = ["AMR_AF_EXAC","AFR_AF_EXAC","NFE_AF_EXAC","FIN_AF_EXAC","OTH_AF_EXAC","GLOBAL_AF_EXAC","EAS_AF_EXAC","SAS_AF_EXAC"]\n')
      fh.write('ops=["concat","concat","concat","concat","concat","concat","concat","concat"]\n\n')
   
   if dbsource == 'oneKG':
      fh.write('fields = ["EAS_AF_1KG","EUR_AF_1KG","AMR_AF_1KG","AFR_AF_1KG","SAS_AF_1KG","GLOBAL_AF_1KG"]\n')
      fh.write('ops=["concat","concat","concat","concat","concat","concat"]\n\n')

   if dbsource == 'docm':
      fh.write('fields = ["DOCM_DISEASE","DOCM_PMID"]\n')
      fh.write('ops=["concat","concat"]\n\n')
      
   if dbsource == 'clinvar':
      fh.write('fields = ["CLINVAR_MSID","CLINVAR_PMIDS","CLINVAR_SIG","CLINVAR_VARIANT_ORIGIN"]\n')
      fh.write('ops=["concat","concat","concat","concat"]\n\n')
   
   fh.close()
   return

def get_vcf_info_tags(vcffile):
   vcf_reader = cyvcf.Reader(open(vcffile, 'r'))
   info_tags = {}
   for info_tag in sorted(vcf_reader.infos.keys()):
      info_tags[str(info_tag)] = 1
   
   return info_tags


def print_vcf_header(query_vcf, vcfheader_file, chromline_only = False):
   if chromline_only == True:
      os.system('bgzip -dc ' + str(query_vcf) + ' | egrep \'^#\' | egrep \'^#CHROM\' >> ' + str(vcfheader_file))
   else:
      os.system('bgzip -dc ' + str(query_vcf) + ' | egrep \'^#\' | egrep -v \'^#CHROM\' > ' + str(vcfheader_file))

if __name__=="__main__": __main__()
