#!/usr/bin/env python

import csv
import re
import argparse
from cyvcf2 import VCF, Writer
import gzip
import os
import annoutils
import oncogenicity

logger = annoutils.getlogger('gvanno-gene-annotate')
csv.field_size_limit(500 * 1024 * 1024)


def __main__():
   
   parser = argparse.ArgumentParser(description='Gene annotations from gvanno pipeline (SNVs/InDels)')
   parser.add_argument('vcf_file', help='VCF file with VEP-annotated query variants (SNVs/InDels)')
   parser.add_argument('gvanno_db_dir',help='gvanno data directory')
   parser.add_argument('lof_prediction',default=0,type=int,help='VEP LoF prediction setting (0/1)')
   parser.add_argument('oncogenicity_annotation',default=0,type=int,help='Include oncogenicity annotation (0/1)')
   parser.add_argument('regulatory_annotation',default=0,type=int,help='Inclusion of VEP regulatory annotations (0/1)')
   args = parser.parse_args()

   extend_vcf_annotations(args.vcf_file, args.gvanno_db_dir, args.lof_prediction, args.oncogenicity_annotation, args.regulatory_annotation)

def extend_vcf_annotations(query_vcf, gvanno_db_directory, lof_prediction = 0, oncogenicity_annotation = 0, regulatory_annotation = 0):
   """
   Function that reads VEP/vcfanno-annotated VCF and extends the VCF INFO column with tags from
   1. CSQ elements within the primary transcript consequence picked by VEP, e.g. SYMBOL, Feature, Gene, Consequence etc.
   2. Gene annotations, e.g. known oncogenes/tumor suppressors
   3. Cancer hotspot mutation sites
   4. Variant effect predictions
   """

   ## read VEP and PCGR tags to be appended to VCF file
   vcf_infotags_meta = annoutils.read_infotag_file(logger, os.path.join(gvanno_db_directory,'gvanno_infotags.tsv'))
   gvanno_xref_map = annoutils.read_genexref_namemap(logger, os.path.join(gvanno_db_directory,'gvanno_xref', 'gvanno_xref_namemap.tsv'))
   out_vcf = re.sub(r'\.vcf(\.gz){0,}$','.annotated.vcf',query_vcf)

   meta_vep_dbnsfp_info = annoutils.vep_dbnsfp_meta_vcf(query_vcf, vcf_infotags_meta)
   dbnsfp_prediction_algorithms = meta_vep_dbnsfp_info['dbnsfp_prediction_algorithms']
   vep_csq_fields_map = meta_vep_dbnsfp_info['vep_csq_fieldmap']
   vcf = VCF(query_vcf)
   for tag in vcf_infotags_meta:
      if lof_prediction == 0 and regulatory_annotation == 0:
         if not tag.startswith('LoF') and not tag.startswith('REGULATORY_'):
            vcf.add_info_to_header({'ID': tag, 'Description': str(vcf_infotags_meta[tag]['description']),'Type':str(vcf_infotags_meta[tag]['type']), 'Number': str(vcf_infotags_meta[tag]['number'])})
      elif lof_prediction == 1 and regulatory_annotation == 0:
         if not tag.startswith('REGULATORY_'):
            vcf.add_info_to_header({'ID': tag, 'Description': str(vcf_infotags_meta[tag]['description']),'Type':str(vcf_infotags_meta[tag]['type']), 'Number': str(vcf_infotags_meta[tag]['number'])})
      elif lof_prediction == 0 and regulatory_annotation == 1:
         if not tag.startswith('LoF'):
            vcf.add_info_to_header({'ID': tag, 'Description': str(vcf_infotags_meta[tag]['description']),'Type':str(vcf_infotags_meta[tag]['type']), 'Number': str(vcf_infotags_meta[tag]['number'])})
      else:
         vcf.add_info_to_header({'ID': tag, 'Description': str(vcf_infotags_meta[tag]['description']),'Type':str(vcf_infotags_meta[tag]['type']), 'Number': str(vcf_infotags_meta[tag]['number'])})


   w = Writer(out_vcf, vcf)
   current_chrom = None
   num_chromosome_records_processed = 0
   num_records_filtered = 0

   cancer_hotspots = annoutils.read_cancer_hotspots(logger, os.path.join(gvanno_db_directory,'cancer_hotspots', 'cancer_hotspots.tsv'))
   
   vcf_info_element_types = {}
   for e in vcf.header_iter():
      header_element = e.info()
      if 'ID' in header_element and 'HeaderType' in header_element and 'Type' in header_element:
         identifier = str(header_element['ID'])
         fieldtype = str(header_element['Type'])
         vcf_info_element_types[identifier] = fieldtype

   for rec in vcf:
      if current_chrom is None:
         current_chrom = str(rec.CHROM)
         num_chromosome_records_processed = 0
      else:
         if str(rec.CHROM) != current_chrom:
            logger.info('Completed summary of functional annotations for ' + str(num_chromosome_records_processed) + ' variants on chromosome ' + str(current_chrom))
            current_chrom = str(rec.CHROM)
            num_chromosome_records_processed = 0
      if rec.INFO.get('CSQ') is None:
         num_records_filtered = num_records_filtered + 1
         #logger.warning('Variant record ' + str(variant_id) + ' does not have CSQ tag from Variant Effect Predictor (--vep_skip_intergenic or --vep_coding_only turned ON?)  - variant will be skipped')
         continue
      num_chromosome_records_processed += 1
      gvanno_xref = annoutils.make_transcript_xref_map(rec, gvanno_xref_map, xref_tag = "GVANNO_XREF")

      if regulatory_annotation == 1:
         csq_record_results_all = annoutils.parse_vep_csq(rec, gvanno_xref, vep_csq_fields_map, logger, pick_only = False, csq_identifier = 'CSQ')

         if 'picked_gene_csq' in csq_record_results_all:
            vep_csq_records_all = csq_record_results_all['picked_gene_csq']
            rec.INFO['REGULATORY_ANNOTATION'] = annoutils.map_regulatory_variant_annotations(vep_csq_records_all)

      vep_csq_record_results = annoutils.parse_vep_csq(rec, gvanno_xref, vep_csq_fields_map, logger, pick_only = True, csq_identifier = 'CSQ')

      
      principal_hgvsp = '.'
      principal_hgvsc = '.'
      if 'picked_csq' in vep_csq_record_results:
         csq_record = vep_csq_record_results['picked_csq']
         for k in csq_record:
            if k in vcf_info_element_types:
               if vcf_info_element_types[k] == "Flag":
                  #rec.INFO[k] = False
                  if csq_record[k] == "1":
                     rec.INFO[k] = True                                
               else:
                  if not csq_record[k] is None:
                     rec.INFO[k] = csq_record[k]

                     if k == 'HGVSp_short':
                        principal_hgvsp = csq_record[k]
                     
                     if k == 'HGVSc':
                        principal_hgvsc = csq_record[k].split(':')[1]
            #else:
            #   print("missing\t" + str(k))

      if 'all_csq' in vep_csq_record_results:
         rec.INFO['VEP_ALL_CSQ'] = ','.join(vep_csq_record_results['all_csq'])
         annoutils.map_cancer_hotspots(vep_csq_record_results['all_csq'], cancer_hotspots, rec, principal_hgvsp, principal_hgvsc)

      if not rec.INFO.get('DBNSFP') is None:
         annoutils.map_variant_effect_predictors(rec, dbnsfp_prediction_algorithms)

      if oncogenicity_annotation == 1:
         oncogenicity.assign_oncogenicity_evidence(rec, tumortype = "Any")
      w.write_record(rec)
   w.close()
   logger.info('Completed summary of functional annotations for ' + str(num_chromosome_records_processed) + ' variants on chromosome ' + str(current_chrom))
   vcf.close()
   logger.info("Number of variant calls filtered by VEP (No CSQ tag, '--vep_coding_only' / '--vep_skip_intergenic'): " + str(num_records_filtered))

   if os.path.exists(out_vcf):
      if os.path.getsize(out_vcf) > 0:
         os.system('bgzip -f ' + str(out_vcf))
         os.system('tabix -f -p vcf ' + str(out_vcf) + '.gz')
         annotated_vcf = out_vcf + '.gz'
         annoutils.write_pass_vcf(annotated_vcf, logger)
      else:
         annoutils.error_message('No remaining PASS variants found in query VCF - exiting and skipping STEP 4 (gvanno-writer)', logger)
   else:
      annoutils.error_message('No remaining PASS variants found in query VCF - exiting and skipping STEP 4 (gvanno-writer)', logger)

if __name__=="__main__": __main__()


      
