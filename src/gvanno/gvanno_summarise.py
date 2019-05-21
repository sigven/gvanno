#!/usr/bin/env python

import csv
import re
import argparse
from cyvcf2 import VCF, Writer
import gzip
import os
import annoutils

logger = annoutils.getlogger('gvanno-gene-annotate')
csv.field_size_limit(500 * 1024 * 1024)


def __main__():
   
   parser = argparse.ArgumentParser(description='Gene annotations from gvanno pipeline (SNVs/InDels)')
   parser.add_argument('vcf_file', help='VCF file with VEP-annotated query variants (SNVs/InDels)')
   parser.add_argument('gvanno_db_dir',help='gvanno data directory')
   parser.add_argument('lof_prediction',default=0,type=int,help='VEP LoF prediction setting (0/1)')
   args = parser.parse_args()

   extend_vcf_annotations(args.vcf_file, args.gvanno_db_dir, args.lof_prediction)

def extend_vcf_annotations(query_vcf, gvanno_db_directory, lof_prediction = 0):
   """
   Function that reads VEP/vcfanno-annotated VCF and extends the VCF INFO column with tags from
   1. CSQ elements within the primary transcript consequence picked by VEP, e.g. SYMBOL, Feature, Gene, Consequence etc.
   2. Gene annotations, e.g. known oncogenes/tumor suppressors, curated disease associations (DisGenet), MIM phenotype associations etc
   3. Protein-relevant annotations, e.g. c functional protein features etc.
   4. Variant effect predictions
   """

   ## read VEP and PCGR tags to be appended to VCF file
   vcf_infotags_meta = annoutils.read_infotag_file(os.path.join(gvanno_db_directory,'gvanno_infotags.tsv'))
   out_vcf = re.sub(r'\.vcf(\.gz){0,}$','.annotated.vcf',query_vcf)

   meta_vep_dbnsfp_info = annoutils.vep_dbnsfp_meta_vcf(query_vcf, vcf_infotags_meta)
   dbnsfp_prediction_algorithms = meta_vep_dbnsfp_info['dbnsfp_prediction_algorithms']
   vep_csq_fields_map = meta_vep_dbnsfp_info['vep_csq_fieldmap']
   vcf = VCF(query_vcf)
   for tag in vcf_infotags_meta:
      if lof_prediction == 0:
         if not tag.startswith('LoF'):
            vcf.add_info_to_header({'ID': tag, 'Description': str(vcf_infotags_meta[tag]['description']),'Type':str(vcf_infotags_meta[tag]['type']), 'Number': str(vcf_infotags_meta[tag]['number'])})
      else:
         vcf.add_info_to_header({'ID': tag, 'Description': str(vcf_infotags_meta[tag]['description']),'Type':str(vcf_infotags_meta[tag]['type']), 'Number': str(vcf_infotags_meta[tag]['number'])})

   
   w = Writer(out_vcf, vcf)
   current_chrom = None
   num_chromosome_records_processed = 0
   gvanno_xref_map = {'ENSEMBL_TRANSCRIPT_ID':0, 'ENSEMBL_GENE_ID':1, 'ENSEMBL_PROTEIN_ID':2, 'SYMBOL':3, 'SYMBOL_ENTREZ':4,'ENTREZ_ID':5, 'UNIPROT_ID':6, 'APPRIS':7,'UNIPROT_ACC':8,
                        'REFSEQ_MRNA':9, 'CORUM_ID':10,'TUMOR_SUPPRESSOR':11,'ONCOGENE':12,'DISGENET_CUI':13,'MIM_PHENOTYPE_ID':14, 'OPENTARGETS_DISEASE_ASSOCS':15,
                        'OPENTARGETS_TRACTABILITY_COMPOUND':16, 'OPENTARGETS_TRACTABILITY_ANTIBODY':17}
   
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
         alt_allele = ','.join(rec.ALT)
         pos = rec.start + 1
         variant_id = 'g.' + str(rec.CHROM) + ':' + str(pos) + str(rec.REF) + '>' + alt_allele
         logger.warning('Variant record ' + str(variant_id) + ' does not have CSQ tag from Variant Effect Predictor (vep_skip_intergenic in config set to true?)  - variant will be skipped')
         continue
      num_chromosome_records_processed += 1
      gvanno_xref = annoutils.make_transcript_xref_map(rec, gvanno_xref_map, xref_tag = "GVANNO_XREF")

      csq_record_results = annoutils.parse_vep_csq(rec, gvanno_xref, vep_csq_fields_map, logger, pick_only = True, csq_identifier = 'CSQ')
      if 'vep_all_csq' in csq_record_results:
         rec.INFO['VEP_ALL_CSQ'] = ','.join(csq_record_results['vep_all_csq'])
      if 'vep_block' in csq_record_results:
         vep_csq_records = csq_record_results['vep_block']
         block_idx = 0
         record = vep_csq_records[block_idx]
         for k in record:
            if k in vcf_info_element_types:
               if vcf_info_element_types[k] == "Flag" and record[k] == "1":
                  rec.INFO[k] = True
               else:
                  if not record[k] is None:
                     rec.INFO[k] = record[k]

      if not rec.INFO.get('DBNSFP') is None:
         annoutils.map_variant_effect_predictors(rec, dbnsfp_prediction_algorithms)

      w.write_record(rec)
   w.close()
   logger.info('Completed summary of functional annotations for ' + str(num_chromosome_records_processed) + ' variants on chromosome ' + str(current_chrom))
   vcf.close()

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


      
