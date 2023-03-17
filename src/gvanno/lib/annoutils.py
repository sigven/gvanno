#!/usr/bin/env python

import os,re,sys
import csv
import logging
import gzip
from cyvcf2 import VCF, Writer
import subprocess


csv.field_size_limit(500 * 1024 * 1024)
threeLettertoOneLetterAA = {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Glu':'E','Gln':'Q','Gly':'G','His':'H',
                           'Ile':'I','Leu':'L','Lys':'K', 'Met':'M','Phe':'F','Pro':'P','Ser':'S','Thr':'T','Trp':'W',
                           'Tyr':'Y','Val':'V','Ter':'X'}



def read_genexref_namemap(logger, gene_xref_namemap_tsv):
   namemap_xref = {} ##dictionary returned
   if not os.path.exists(gene_xref_namemap_tsv):
      logger.info("ERROR: File '" + str(gene_xref_namemap_tsv) + "' does not exist - exiting")
      exit(1)
      #return namemap_xref
   tsvfile = open(gene_xref_namemap_tsv, 'r')
   reader = csv.DictReader(tsvfile, delimiter='\t')
   for row in reader:
      namemap_xref[row['name']] = int(row['index'])
   tsvfile.close()

   return namemap_xref

def read_infotag_file(logger, vcf_info_tags_tsv):
   """
   Function that reads a file that lists VCF INFO tags produced by PCGR/CPSR/gvanno.
   An example line of the VCF info tag file is the following:
   
   tag	number	type	description category
   Consequence	.	String	"Impact modifier for the consequence type (picked by VEP's --flag_pick_allele option)."   vep
   
   A dictionary is returned, with the tag as the key, and the full dictionary record as the value
   """
   info_tag_xref = {} ##dictionary returned
   if not os.path.exists(vcf_info_tags_tsv):
      logger.info("ERROR: File '" + str(vcf_info_tags_tsv) + "' does not exist - exiting")
      return info_tag_xref
   tsvfile = open(vcf_info_tags_tsv, 'r')
   reader = csv.DictReader(tsvfile, delimiter='\t')
   for row in reader:
      if not row['tag'] in info_tag_xref:
         info_tag_xref[row['tag']] = row
   tsvfile.close()
   
   return info_tag_xref

def check_subprocess(command):
   # if debug:
   #    logger.info(command)
   try:
      output = subprocess.check_output(str(command), stderr=subprocess.STDOUT, shell=True)
      if len(output) > 0:
         print (str(output.decode()).rstrip())
   except subprocess.CalledProcessError as e:
      print (e.output.decode())
      exit(0)


def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()

def error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
   exit(1)

def write_pass_vcf(annotated_vcf, logger):
   """
   Function prints all PASS variants from a VCF file. New VCF file appends '.pass.' to the filename.
   
   """
   out_vcf = re.sub(r'\.vcf\.gz$','.pass.vcf',annotated_vcf)
   vcf = VCF(annotated_vcf)
   w = Writer(out_vcf, vcf)

   num_rejected = 0
   num_pass = 0
   for rec in vcf:
      if rec.FILTER is None or rec.FILTER == 'None':
         w.write_record(rec)
         num_pass += 1
      else:
         num_rejected +=1

   vcf.close()
   w.close()
   
   logger.info('Number of non-PASS/REJECTED variant calls: ' + str(num_rejected))
   logger.info('Number of PASSed variant calls: ' + str(num_pass))
   if num_pass == 0:
      logger.warning('There are zero variants with a \'PASS\' filter in the VCF file')
      os.system('bgzip -dc ' + str(annotated_vcf) + ' egrep \'^#\' > ' + str(out_vcf))
   #else:
   os.system('bgzip -f ' + str(out_vcf))
   os.system('tabix -f -p vcf ' + str(out_vcf) + '.gz')

   return

def warn_message(message, logger):
   logger.warning(message)

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


def is_valid_vcf(input_vcf, output_dir, logger, debug):
   """
   Function that reads the output file of EBIvariation/vcf-validator and reports potential errors and validation status
   """

   logger.info('Validating VCF file with EBIvariation/vcf-validator (v0.9.3)')

   ## Filename for output file from VCF validation
   vcf_validation_output_file = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)', '.vcf_validator_output', os.path.basename(input_vcf)))
   
   ## Command for running vcf-validator
   command_v42 = 'vcf_validator -i ' + str(input_vcf) + ' -l warning -r text -o ' + str(output_dir) + ' > ' + str(vcf_validation_output_file)  + ' 2>&1'
   if debug:
      logger.info(command_v42)
   try:
      subprocess.run(command_v42, stderr=subprocess.STDOUT, shell=True)
   except subprocess.CalledProcessError as e:
      print (e.output.decode())
      exit(0)

   ## read output file from vcf-validator and get name of log/error file
   validation_results = {}
   validation_results['validation_status'] = 0
   validation_results['error_messages'] = []
   validation_log_fname = None
   if os.path.exists(vcf_validation_output_file):
      f = open(vcf_validation_output_file, 'r')
      for line in f:
         if 'Text report written to' in line.rstrip():
            validation_log_fname = re.split(r' : ',line.rstrip(),maxsplit=2)[1]
      f.close()
      if not debug:
         check_subprocess('rm -f ' + str(vcf_validation_output_file))

   if validation_log_fname is None:
      err_msg = 'Cannot find file with error messages from vcf_validator'
      return error_message(err_msg, logger)

   if os.path.exists(validation_log_fname):
      f = open(validation_log_fname, 'r')
      for line in f:
         if not re.search(r' \(warning\)$|^Reading from ',line.rstrip()): ## ignore warnings
            if line.startswith('Line '):
               validation_results['error_messages'].append('ERROR: ' + line.rstrip())
            if 'the input file is valid' in line.rstrip(): ## valid VCF
               validation_results['validation_status'] = 1
            if 'the input file is not valid' in line.rstrip():  ## non-valid VCF
               validation_results['validation_status'] = 0
      f.close()
      if not debug:
         check_subprocess('rm -f ' + str(validation_log_fname))
   else:
      err_msg = str(validation_log_fname) + ' does not exist'
      return error_message(err_msg, logger)

   if validation_results['validation_status'] == 0:
      error_string_42 = '\n'.join(validation_results['error_messages'])
      validation_status = 'According to the VCF specification, the VCF file (' + str(input_vcf) + ') is NOT valid'
      err_msg = validation_status + ':\n' + str(error_string_42)
      return error_message(err_msg, logger)
   else:
      validation_status = 'According to the VCF specification, the VCF file ' + str(input_vcf) + ' is valid'
      logger.info(validation_status)
   return 0


def map_regulatory_variant_annotations(vep_csq_records):
   """
   Function that considers an array of VEP CSQ records and appends all regulatory variant consequent annotations (open chromatin, TF_binding_site, 
   CTCF_binding_site, promoter (flanks), enhancers ) into a single comma-separated string. Each individual regulatory annotation is formatted as:
   <Consequence>|<Feature_type>|<Feature>|<BIOTYPE>|<MOTIF_NAME>|<MOTIF_POS>|<HIGH_INF_POS>|<MOTIF_SCORE_CHANGE>|<TRANSCRIPTION_FACTORS>
   """

   regulatory_annotation = '.'
   if len(vep_csq_records) == 1:
      return regulatory_annotation

   j = 0
   regulatory_annotations = []
   while j < len(vep_csq_records):
      missing_key_annotation = False
      for k in ['Feature_type','Consequence','Feature']:
         if not k in vep_csq_records[j].keys():
            missing_key_annotation = True

      if missing_key_annotation is False:

         ## RegulatoryFeature annotations - open chromatin, promoters (flanks), enhancers, CTCF binding sites
         if vep_csq_records[j]['Feature_type'] == 'RegulatoryFeature':
            biotype = ""
            if re.match(r"^(enhancer|promoter|open|CTCF|TF_)", vep_csq_records[j]['BIOTYPE']):
               biotype = vep_csq_records[j]['BIOTYPE']
            
            annotation = str(vep_csq_records[j]['Consequence']) + '|' + \
               str(vep_csq_records[j]['Feature_type']) + '|' + \
               str(vep_csq_records[j]['Feature']) + '|' + \
               str(biotype) + '|||||'

            regulatory_annotations.append(annotation)

         ## MotifFeature annotations (TF)
         if vep_csq_records[j]['Feature_type'] == 'MotifFeature':
            missing_motif_annotation = False
            for annotation in ['MOTIF_NAME','MOTIF_POS','HIGH_INF_POS','MOTIF_SCORE_CHANGE','TRANSCRIPTION_FACTORS']:
               if not annotation in vep_csq_records[j].keys():
                  missing_motif_annotation = True
            
            if missing_motif_annotation is False:
               annotation = str(vep_csq_records[j]['Consequence']) + '|' + \
                  str(vep_csq_records[j]['Feature_type']) + '|' + \
                  str(vep_csq_records[j]['Feature']) + '|TF_binding_site|' + \
                  str(vep_csq_records[j]['MOTIF_NAME']) + '|' + \
                  str(vep_csq_records[j]['MOTIF_POS']) + '|' + \
                  str(vep_csq_records[j]['HIGH_INF_POS']) + '|' + \
                  str(vep_csq_records[j]['MOTIF_SCORE_CHANGE']) + '|' + \
                  str(vep_csq_records[j]['TRANSCRIPTION_FACTORS'])

               regulatory_annotations.append(annotation)
      
      j = j + 1
   
   if len(regulatory_annotations) > 0:
      regulatory_annotation = ','.join(regulatory_annotations)

   return(regulatory_annotation)
               
def get_correct_cpg_transcript(vep_csq_records):
   """
   Function that considers an array of VEP CSQ records and picks most relevant consequence (and gene) from
   neighbouring genes/transcripts of relevance for cancer predisposition (cpg = cancer predisposition gene)
   """

  
   csq_idx = 0
   if len(vep_csq_records) == 1:
      return csq_idx
   
   
   ## some variants are assigned multiple transcript consequences
   ## if cancer predisposition genes are in the vicinity of other genes, choose the cancer predisposition gene
   ## if there are neighbouring cancer-predispositon genes, choose custom gene, preferring coding change (see below, KLLN/PTEN, XPC/TMEM43, NTHL1/TSC2)
   csq_idx_dict = {}
   for g in ['KLLN','PTEN','XPC','TMEM43','NTHL1','TSC2']:
      csq_idx_dict[g] = {}
      csq_idx_dict[g]['idx'] = -1
      csq_idx_dict[g]['coding'] = False

   j = 0
   while j < len(vep_csq_records):
      if 'CANCER_PREDISPOSITION_SOURCE' in vep_csq_records[j].keys() or 'GE_PANEL_ID' in vep_csq_records[j].keys():
         csq_idx = j
         if 'SYMBOL' in vep_csq_records[j].keys():
            if vep_csq_records[j]['SYMBOL'] in csq_idx_dict.keys():
               csq_idx_dict[str(vep_csq_records[j]['SYMBOL'])]['idx'] = j
               if vep_csq_records[j]['CODING_STATUS'] == 'coding':
                  csq_idx = j  # prefer coding on over anything else
                  csq_idx_dict[str(vep_csq_records[j]['SYMBOL'])]['coding'] = True
      j = j + 1
   
   if csq_idx_dict['KLLN']['idx'] != -1 and csq_idx_dict['PTEN']['idx'] != -1:
      csq_idx = csq_idx_dict['PTEN']['idx']
      if csq_idx_dict['KLLN']['coding'] is True:
         csq_idx = csq_idx_dict['KLLN']['idx']
   
   if csq_idx_dict['XPC']['idx'] != -1 and csq_idx_dict['TMEM43']['idx'] != -1:
      csq_idx = csq_idx_dict['XPC']['idx']
      if csq_idx_dict['TMEM43']['coding'] is True:
         csq_idx = csq_idx_dict['TMEM43']['idx']
   
   if csq_idx_dict['TSC2']['idx'] != -1 and csq_idx_dict['NTHL1']['idx'] != -1:
      csq_idx = csq_idx_dict['TSC2']['idx']
      if csq_idx_dict['NTHL1']['coding'] is True:
         csq_idx = csq_idx_dict['NTHL1']['idx']

   if csq_idx is None:
      csq_idx = 0
   return csq_idx


def threeToOneAA(aa_change):
	
   for three_letter_aa in threeLettertoOneLetterAA.keys():
      aa_change = aa_change.replace(three_letter_aa,threeLettertoOneLetterAA[three_letter_aa])

   aa_change = re.sub(r'[A-Z]{1}fsX([0-9]{1,}|\?)','fs',aa_change)
   return aa_change

def map_variant_effect_predictors(rec, algorithms):
   
   dbnsfp_predictions = map_dbnsfp_predictions(str(rec.INFO.get('DBNSFP')), algorithms)
   if rec.INFO.get('Gene') is None or rec.INFO.get('Consequence') is None:
      return
   gene_id = str(rec.INFO.get('Gene'))
   consequence = str(rec.INFO.get('Consequence'))
     
   dbnsfp_key = ''

   found_key = 0
   if not rec.INFO.get('HGVSp_short') is None and not rec.INFO.get('HGVSp_short') == '.':
      aa_change = str(rec.INFO.get('HGVSp_short'))
      dbnsfp_key = gene_id + ':' + str(aa_change)
      if dbnsfp_key in dbnsfp_predictions:
         found_key = 1
   
   if found_key == 0 and re.search('splice_',consequence):
      dbnsfp_key = gene_id

   if dbnsfp_key != '':
      if dbnsfp_key in dbnsfp_predictions:
         rec.INFO['EFFECT_PREDICTIONS'] = dbnsfp_predictions[dbnsfp_key]
         for algo_pred in rec.INFO['EFFECT_PREDICTIONS'].split('&'):
            if algo_pred.startswith('sift:'):
               rec.INFO['DBNSFP_SIFT'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('provean:'):
               rec.INFO['DBNSFP_PROVEAN'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('m-cap:'):
               rec.INFO['DBNSFP_M_CAP'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('mutpred:'):
               rec.INFO['DBNSFP_MUTPRED'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('metarnn:'):
               rec.INFO['DBNSFP_META_RNN'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('fathmm:'):
               rec.INFO['DBNSFP_FATHMM'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('fathmm_mkl_coding:'):
               rec.INFO['DBNSFP_FATHMM_MKL'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('mutationtaster:'):
               rec.INFO['DBNSFP_MUTATIONTASTER'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('mutationassessor:'):
               rec.INFO['DBNSFP_MUTATIONASSESSOR'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('deogen2:'):
               rec.INFO['DBNSFP_DEOGEN2'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('primateai:'):
               rec.INFO['DBNSFP_PRIMATEAI'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('list_s2:'):
               rec.INFO['DBNSFP_LIST_S2'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('gerp_rs:'):
               rec.INFO['DBNSFP_GERP'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('bayesdel_addaf:'):
               rec.INFO['DBNSFP_BAYESDEL_ADDAF'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('aloft:'):
               rec.INFO['DBNSFP_ALOFTPRED'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('splice_site_rf:'):
               rec.INFO['DBNSFP_SPLICE_SITE_RF'] = str(algo_pred.split(':')[1])
            if algo_pred.startswith('splice_site_ada:'):
               rec.INFO['DBNSFP_SPLICE_SITE_ADA'] = str(algo_pred.split(':')[1])



def detect_reserved_info_tag(tag, tag_name, logger):
   reserved_tags = ['AA','AC','AF','AN','BQ','CIGAR','DB','DP','END','H2','H3','MQ','MQ0','NS','SB','SOMATIC','VALIDATED','1000G']
   if tag in reserved_tags:
      err_msg = 'Custom INFO tag (' + str(tag_name) + ') needs another name - ' + str(tag) + ' is a reserved field in the VCF specification (INFO)'
      return error_message(err_msg, logger)
   
   reserved_format_tags = ['GT','DP','FT','GL','GLE','GQ','PL','HQ','PS','PQ','EC','MQ']
   if tag in reserved_format_tags:
      err_msg = 'Custom INFO tag (' + str(tag_name) + ') needs another name - ' + str(tag) + ' is a reserved field in the VCF specification (FORMAT)'
      return error_message(err_msg, logger)


def assign_cds_exon_intron_annotations(csq_record):
   csq_record['CODING_STATUS'] = 'noncoding'
   csq_record['EXONIC_STATUS'] = 'nonexonic'
   csq_record['SPLICE_DONOR_RELEVANT'] = False
   csq_record['NULL_VARIANT'] = False
   csq_record['INTRON_POSITION'] = 0
   csq_record['EXON_POSITION'] = 0
   
   coding_csq_pattern = r"^(stop_|start_lost|frameshift_|missense_|splice_donor|splice_acceptor|protein_altering|inframe_)"
   wes_csq_pattern = r"^(stop_|start_lost|frameshift_|missense_|splice_donor|splice_acceptor|inframe_|protein_altering|synonymous)"
   null_pattern = r"^(stop_|frameshift_)"
   if re.match(coding_csq_pattern, str(csq_record['Consequence'])):
      csq_record['CODING_STATUS'] = 'coding'
   
   if re.match(wes_csq_pattern, str(csq_record['Consequence'])):
      csq_record['EXONIC_STATUS'] = 'exonic'
   
   if re.match(null_pattern, str(csq_record['Consequence'])):
      csq_record['NULL_VARIANT'] = True

   if re.match(r"^splice_region_variant", str(csq_record['Consequence'])) and re.search(r'(\+3(A|G)>|\+4G>|\+5G>)', str(csq_record['HGVSc'])):
      csq_record['SPLICE_DONOR_RELEVANT'] = True

   if re.match(r"splice_region_variant|intron_variant", str(csq_record['Consequence'])):
      match = re.search(r"((-|\+)[0-9]{1,}(dup|del|inv|((ins|del|dup|inv|delins)(A|G|C|T){1,})|(A|C|T|G){1,}>(A|G|C|T){1,}))$", str(csq_record['HGVSc']))
      if match is not None:
         pos = re.sub(r"(\+|dup|del|delins|ins|inv|(A|G|C|T){1,}|>)","",match.group(0))
         if is_integer(pos):
            csq_record['INTRON_POSITION'] = int(pos)
   
   if 'NearestExonJB' in csq_record.keys():
      if not csq_record['NearestExonJB'] is None:
         if re.match(r"synonymous_|missense_|stop_|inframe_|start_", str(csq_record['Consequence'])) and str(csq_record['NearestExonJB']) != "":
            exon_pos_info = csq_record['NearestExonJB'].split("+")
            if len(exon_pos_info) == 4:
               if is_integer(exon_pos_info[1]) and str(exon_pos_info[2]) == "end":
                  csq_record['EXON_POSITION'] = -int(exon_pos_info[1])
               if is_integer(exon_pos_info[1]) and str(exon_pos_info[2]) == "start":
                  csq_record['EXON_POSITION'] = int(exon_pos_info[1])


   for m in ['HGVSp_short','CDS_CHANGE']:
      csq_record[m] = '.'
   if not csq_record['HGVSc'] is None:
      if csq_record['HGVSc'] != '.':
         if 'splice_acceptor_variant' in csq_record['Consequence'] or 'splice_donor_variant' in csq_record['Consequence']:
            key = str(csq_record['Consequence']) + ':' + str(csq_record['HGVSc'])
            csq_record['CDS_CHANGE'] = key
   if csq_record['Amino_acids'] is None or csq_record['Protein_position'] is None or csq_record['Consequence'] is None:
      return
   if not csq_record['Protein_position'] is None:
      if csq_record['Protein_position'].startswith('-'):
         return

   protein_change = '.'
   if '/' in csq_record['Protein_position']:
      protein_position = str(csq_record['Protein_position'].split('/')[0])
      if '-' in protein_position:
         if protein_position.split('-')[0].isdigit():
            csq_record['AMINO_ACID_START'] = protein_position.split('-')[0]
         if protein_position.split('-')[1].isdigit():
            csq_record['AMINO_ACID_END'] = protein_position.split('-')[1]
      else:
         if protein_position.isdigit():
            csq_record['AMINO_ACID_START'] = protein_position
            csq_record['AMINO_ACID_END'] = protein_position
   
   if not csq_record['HGVSp'] is None:
      if csq_record['HGVSp'] != '.':
         if ':' in csq_record['HGVSp']:
            protein_identifier = str(csq_record['HGVSp'].split(':')[0])
            if protein_identifier.startswith('ENSP'):
               protein_change_VEP = str(csq_record['HGVSp'].split(':')[1])
               protein_change = threeToOneAA(protein_change_VEP)
  
   if 'synonymous_variant' in csq_record['Consequence']:
      protein_change = 'p.' + str(csq_record['Amino_acids']) + str(protein_position) + str(csq_record['Amino_acids'])
      if 'stop_lost' in str(csq_record['Consequence']) and '/' in str(csq_record['Amino_acids']):
         protein_change = 'p.X' + str(protein_position) + str(csq_record['Amino_acids']).split('/')[1]
    
   csq_record['HGVSp_short'] = protein_change
   exon_number = 'NA'
   if not csq_record['EXON'] is None:
      if csq_record['EXON'] != '.':
         if '/' in csq_record['EXON']:
            exon_number = str(csq_record['EXON']).split('/')[0]
            num_exons = str(csq_record['EXON']).split('/')[1]
            if exon_number == num_exons:
               csq_record['LAST_EXON'] = True

   if not csq_record['INTRON'] is None:
      if csq_record['INTRON'] != '.':
         if '/' in csq_record['INTRON']:
            intron_number = str(csq_record['INTRON']).split('/')[0]
            num_introns = str(csq_record['INTRON']).split('/')[1]
            if intron_number == num_introns:
               csq_record['LAST_INTRON'] = True

   if not csq_record['HGVSc'] is None:
      if csq_record['HGVSc'] != '.':
         if protein_change != '.':
            key = str(csq_record['Consequence']) + ':' + str(csq_record['HGVSc']) + ':exon' + str(exon_number) + ':' + str(protein_change)
            csq_record['CDS_CHANGE'] = key

   return

def map_dbnsfp_predictions(dbnsfp_tag, algorithms):
   
   effect_predictions = {}
   
   for v in dbnsfp_tag.split(','):
   
      dbnsfp_info = v.split('|')
      if len(dbnsfp_info) == 1:
         return effect_predictions
      ref_aa = dbnsfp_info[0]
      alt_aa = dbnsfp_info[1]
      all_ids = dbnsfp_info[4].split('&')
      unique_ids = {}
      for s in all_ids:
         unique_ids[s] = 1
         
      isoform_aa_keys = []
      if ref_aa != '.' and alt_aa != '.' and ref_aa != '' and alt_aa != '':
         aa_pos = dbnsfp_info[6].split('&')
         for pos in aa_pos:
            for gene_id in unique_ids:
               k = str(gene_id) + ':p.' + str(ref_aa) + pos + str(alt_aa)
               isoform_aa_keys.append(k)
      else:
         #continue
         for gene_id in unique_ids:
            isoform_aa_keys.append(gene_id)
   
      algorithm_raw_predictions = {}
   
      i = 7
      v = 0
      
      if len(algorithms) != len(dbnsfp_info[7:]):
         return effect_predictions
      
      while i < len(dbnsfp_info):
         algorithm_raw_predictions[str(algorithms[v]).lower()] = dbnsfp_info[i].split('&')
         i = i + 1
         v = v + 1
      dbnsfp_predictions = {}
      
      for k in isoform_aa_keys:
         if not k in dbnsfp_predictions:
            dbnsfp_predictions[k] = {}
         all_preds = []
         for algo in algorithm_raw_predictions.keys():
            unique_algo_predictions = {}
            for pred in algorithm_raw_predictions[algo]:
               if pred != '':
                  if not pred in unique_algo_predictions:
                     unique_algo_predictions[pred] = 1
               else:
                  unique_algo_predictions['.'] = 1
            if len(unique_algo_predictions.keys()) > 1 and '.' in unique_algo_predictions.keys():
               del unique_algo_predictions['.']
            dbnsfp_predictions[k][algo] = str(algo) + ':' + '|'.join(unique_algo_predictions.keys())  
            all_preds.append(dbnsfp_predictions[k][algo])
         effect_predictions[k] = '&'.join(all_preds)

   return effect_predictions


def make_transcript_xref_map(rec, fieldmap, xref_tag = 'PCGR_ONCO_XREF'):
   transcript_xref_map = {}
   if not rec.INFO.get(xref_tag) is None:
      for tref in rec.INFO.get(xref_tag).split(','):
         xrefs = tref.split('|')
         ensembl_transcript_id = str(xrefs[0])
         transcript_xref_map[ensembl_transcript_id] = {}
         for annotation in fieldmap.keys():
            annotation_index = fieldmap[annotation]
            if annotation_index > (len(xrefs) - 1):
               continue
            if xrefs[annotation_index] != '':
               transcript_xref_map[ensembl_transcript_id][annotation] = xrefs[annotation_index]
   
   return(transcript_xref_map)

def vep_dbnsfp_meta_vcf(query_vcf, info_tags_wanted):
   vep_to_pcgr_af = {'gnomAD_AMR_AF':'AMR_AF_GNOMAD',
                     'gnomAD_AFR_AF':'AFR_AF_GNOMAD',
                     'gnomAD_EAS_AF':'EAS_AF_GNOMAD',
                     'gnomAD_NFE_AF':'NFE_AF_GNOMAD',
                     'gnomAD_AF':'GLOBAL_AF_GNOMAD',
                     'gnomAD_SAS_AF':'SAS_AF_GNOMAD',
                     'gnomAD_OTH_AF':'OTH_AF_GNOMAD',
                     'gnomAD_ASJ_AF':'ASJ_AF_GNOMAD',
                     'gnomAD_FIN_AF':'FIN_AF_GNOMAD',
                     'AFR_AF':'AFR_AF_1KG',
                     'AMR_AF':'AMR_AF_1KG',
                     'SAS_AF':'SAS_AF_1KG',
                     'EUR_AF':'EUR_AF_1KG',
                     'EAS_AF':'EAS_AF_1KG', 
                     'AF':'GLOBAL_AF_1KG'}

   vcf = VCF(query_vcf)
   vep_csq_index2fields = {}
   vep_csq_fields2index = {}
   dbnsfp_prediction_algorithms = []
   for e in vcf.header_iter():
      header_element = e.info()
      if 'ID' in header_element.keys():
         identifier = str(header_element['ID'])
         if identifier == 'CSQ' or identifier == 'DBNSFP':
            description = str(header_element['Description'])
            if 'Format: ' in description:
               subtags = description.split('Format: ')[1].split('|')
               if identifier == 'CSQ':
                  i = 0
                  for t in subtags:
                     v = t.replace('"','')
                     if t in vep_to_pcgr_af:
                        v = str(vep_to_pcgr_af[t])
                     if v in info_tags_wanted:
                        vep_csq_index2fields[i] = v
                        vep_csq_fields2index[v] = i
                     i = i + 1
               if identifier == 'DBNSFP':
                  i = 7
                  while(i < len(subtags)):
                     dbnsfp_prediction_algorithms.append(str(re.sub(r'((_score)|(_pred))"*$','',subtags[i])))
                     i = i + 1

   vep_dbnsfp_meta_info = {}
   vep_dbnsfp_meta_info['vep_csq_fieldmap'] = {}
   vep_dbnsfp_meta_info['vep_csq_fieldmap']['field2index'] = vep_csq_fields2index
   vep_dbnsfp_meta_info['vep_csq_fieldmap']['index2field'] = vep_csq_index2fields
   vep_dbnsfp_meta_info['dbnsfp_prediction_algorithms'] = dbnsfp_prediction_algorithms

   return vep_dbnsfp_meta_info

def parse_vep_csq(rec, transcript_xref_map, vep_csq_fields_map, logger, pick_only = True, csq_identifier = 'CSQ'):

   all_csq_pick = []
   all_transcript_consequences = []
   
   for csq in rec.INFO.get(csq_identifier).split(','):
      csq_fields =  csq.split('|')

      entrezgene = '.'

      ## Entrez gene identifier is not provided directly by VEP, pull out this from 'transcript_xref_map' for a given transcript-specific CSQ block
      ##  - used for 'consequence_entry' object that are added to 'vep_all_csq' array
      k = 0
      while(k < len(csq_fields)):
         if k in vep_csq_fields_map['index2field']:
            if vep_csq_fields_map['index2field'][k] == 'Feature':
               ensembl_transcript_id = csq_fields[k]
               if ensembl_transcript_id != '' and ensembl_transcript_id.startswith('ENST'):
                  if ensembl_transcript_id in transcript_xref_map.keys():
                     if 'ENTREZGENE' in transcript_xref_map[ensembl_transcript_id].keys():
                        entrezgene = transcript_xref_map[ensembl_transcript_id]['ENTREZGENE']
         k = k + 1


      if pick_only is False:
         j = 0
         csq_record = {}

         while(j < len(csq_fields)):
            if j in vep_csq_fields_map['index2field']:
               if csq_fields[j] != '':
                  csq_record[vep_csq_fields_map['index2field'][j]] = str(csq_fields[j])
            j = j + 1
         all_csq_pick.append(csq_record)

      else:
         ## loop over VEP consequence blocks PICK'ed according to VEP's ranking scheme
         if csq_fields[vep_csq_fields_map['field2index']['PICK']] == "1": ## only consider the primary/picked consequence when expanding with annotation tags
            j = 0
            csq_record = {}

            ## loop over block annotation elements (separated with '|'), and assign their values to the csq_record dictionary object
            while(j < len(csq_fields)):
               if j in vep_csq_fields_map['index2field']:
                  if csq_fields[j] != '':
                     csq_record[vep_csq_fields_map['index2field'][j]] = str(csq_fields[j])
                     if vep_csq_fields_map['index2field'][j] == 'Feature':
                        ensembl_transcript_id = str(csq_fields[j])
                        if ensembl_transcript_id in transcript_xref_map:
                           for annotation in transcript_xref_map[ensembl_transcript_id].keys():
                              if annotation != 'SYMBOL':
                                 ## assign additional gene/transcript annotations from the custom transcript 
                                 ## xref map (PCGR/CPSR/gvanno) as key,value pairs in the csq_record object
                                 csq_record[annotation] = transcript_xref_map[ensembl_transcript_id][annotation]
                                 #print(str(annotation) + "\t" + str(transcript_xref_map[ensembl_transcript_id][annotation]))
                                 
                              else:
                                 #print(str(annotation) + "\t" + str(transcript_xref_map[ensembl_transcript_id][annotation]))
                                 ## If not VEP provides SYMBOL, append SYMBOL provided by xref map
                                 if csq_record['SYMBOL'] is None:
                                     csq_record[annotation] = transcript_xref_map[ensembl_transcript_id][annotation]
                        else:
                           if re.match(r'ENST', ensembl_transcript_id):
                              logger.warning('Could not find transcript xrefs for ' + str(ensembl_transcript_id))

                     ## Specifically assign PFAM protein domain as a csq_record key
                     if vep_csq_fields_map['index2field'][j] == 'DOMAINS':
                        domain_identifiers = str(csq_fields[j]).split('&')
                        for v in domain_identifiers:
                           if v.startswith('Pfam'):
                              csq_record['PFAM_DOMAIN'] = str(re.sub(r'\.[0-9]{1,}$','',re.sub(r'Pfam:','',v)))
                        
                        csq_record['DOMAINS'] = None
                     ## Assign COSMIC/DBSNP mutation ID's as individual key,value pairs in the csq_record object
                     if vep_csq_fields_map['index2field'][j] == 'Existing_variation':
                        var_identifiers = str(csq_fields[j]).split('&')
                        cosmic_identifiers = []
                        dbsnp_identifiers = []
                        for v in var_identifiers:
                           if v.startswith('COSV'):
                              cosmic_identifiers.append(v)
                           if v.startswith('COSM'):
                              cosmic_identifiers.append(v)
                           if v.startswith('rs'):
                              dbsnp_identifiers.append(v)
                        if len(cosmic_identifiers) > 0:
                           csq_record['COSMIC_MUTATION_ID'] = '&'.join(cosmic_identifiers)
                        if len(dbsnp_identifiers) > 0:
                           csq_record['DBSNPRSID'] = '&'.join(dbsnp_identifiers)
                  else:
                     csq_record[vep_csq_fields_map['index2field'][j]] = None
               j = j + 1
            
            ## Assign coding status, protein change, coding sequence change, last exon/intron status etc
            assign_cds_exon_intron_annotations(csq_record)
            ## Append transcript consequence to all_csq_pick
            all_csq_pick.append(csq_record)
         symbol = "."
         hgvsc = "."
         hgvsp = "."
         if csq_fields[vep_csq_fields_map['field2index']['SYMBOL']] != "":
            symbol = str(csq_fields[vep_csq_fields_map['field2index']['SYMBOL']])
         if csq_fields[vep_csq_fields_map['field2index']['HGVSc']] != "":
            hgvsc = str(csq_fields[vep_csq_fields_map['field2index']['HGVSc']].split(':')[1])
         if csq_fields[vep_csq_fields_map['field2index']['HGVSp']] != "":
            hgvsp = str(csq_fields[vep_csq_fields_map['field2index']['HGVSp']].split(':')[1])
         consequence_entry = (str(csq_fields[vep_csq_fields_map['field2index']['Consequence']]) + ":" +  
            str(symbol) + ":" + 
            str(entrezgene) + ":" +
            str(hgvsc) + ":" + 
            str(hgvsp) + ":" + 
            str(csq_fields[vep_csq_fields_map['field2index']['Feature_type']]) + ":" + 
            str(csq_fields[vep_csq_fields_map['field2index']['Feature']]) + ":" + 
            str(csq_fields[vep_csq_fields_map['field2index']['BIOTYPE']]))
         all_transcript_consequences.append(consequence_entry)


   vep_chosen_csq_idx = 0
   vep_csq_results = {}
   vep_csq_results['picked_gene_csq'] = all_csq_pick
   vep_csq_results['all_csq'] = all_transcript_consequences
   vep_csq_results['picked_csq'] = None

   ## IF multiple transcript-specific variant consequences highlighted by --pick_allele_gene , 
   ## prioritize/choose block of consequence which has
   ## - A gene with BIOTYPE equal to 'protein-coding' (the other picked transcript/gene may potentialy carry another BIOTYPE nature)
   ## - A gene consequence classified as 'exonic' (the other picked transcript/gene likely carries a nonexonic consequence)
   if len(vep_csq_results['picked_gene_csq']) > 0:
      vep_selected_idx = {}
      vep_selected_idx['exonic_status'] = {}
      vep_selected_idx['consequence'] = {}

      i = 0
      #print('')
      for trans_rec in vep_csq_results['picked_gene_csq']:
         if 'BIOTYPE' in trans_rec and 'Consequence' in trans_rec and 'EXONIC_STATUS' in trans_rec:
            if not trans_rec['BIOTYPE'] is None and not trans_rec['Consequence'] is None:
               if trans_rec['BIOTYPE'] == "protein_coding":
                  ## for protein-coding genes - record the exonic variant consequence status (exonic/nonexonic)
                  vep_selected_idx['exonic_status'][i] = trans_rec['EXONIC_STATUS']
                  vep_selected_idx['consequence'][i] = trans_rec['Consequence']
         i = i + 1
            
      ## when multiple transcript gene blocks are picked by VEP, prioritize the block with 'exonic' consequence
      if len(vep_selected_idx['exonic_status'].keys()) > 1:
         exonic_cons_found = 0
         for j in vep_selected_idx['exonic_status'].keys():
            if vep_selected_idx['exonic_status'][j] == 'exonic':
               exonic_cons_found = 1
               vep_chosen_csq_idx = j

         ## if multiple non-exonic variants are found, prioritize UTR variants over other nonexonic
         if exonic_cons_found == 0:
            for j in vep_selected_idx['consequence'].keys():
               if vep_selected_idx['consequence'][j] == '5_prime_UTR_variant' or vep_selected_idx['consequence'][j] == '3_prime_UTR_variant':                  
                  vep_chosen_csq_idx = j

      else:
         if len(vep_selected_idx['exonic_status'].keys()) == 1:
            for k in vep_selected_idx['exonic_status'].keys():
               vep_chosen_csq_idx = k
      
      vep_csq_results['picked_csq'] = vep_csq_results['picked_gene_csq'][vep_chosen_csq_idx]

   else:
      print('ERROR: No VEP block chosen by --pick_allele_gene')

   return(vep_csq_results)


def read_cancer_hotspots(logger, hotspots_tsv):

   ## load cancer hotspot entries from 'cancer_hotspots.tsv' (provided by github.com/sigven/cancerHotspots)

   hotspots = {} ##dictionary returned
   mutation_hotspots = {}
   codon_hotspots = {}
   splice_hotspots = {}
   if not os.path.exists(hotspots_tsv):
      logger.info("ERROR: File '" + str(hotspots_tsv) + "' does not exist - exiting")
      exit(1)
   tsvfile = open(hotspots_tsv, 'r')
   reader = csv.DictReader(tsvfile, delimiter='\t')
   for row in reader:
      mutation_hotspots[str(row['entrezgene']) + '-' + row['hgvsp2']] = row
      codon_hotspots[str(row['entrezgene']) + '-' + row['codon']] = row
      if row['hgvsc'] != '.':
         splice_hotspots[str(row['entrezgene']) + '-' + str(row['hgvsc'])] = row
   tsvfile.close()

   hotspots['mutation'] = mutation_hotspots
   hotspots['codon'] = codon_hotspots
   hotspots['splice'] = splice_hotspots
   return hotspots


def map_cancer_hotspots(transcript_csq_elements, cancer_hotspots, rec, principal_hgvsp, principal_hgvsc):

   unique_hotspot_mutations = {}
   unique_hotspot_codons = {}

   principal_codon = '.'
   if re.match(r'^(p.[A-Z]{1}[0-9]{1,}[A-Za-z]{1,})', principal_hgvsp):
      codon_match = re.findall(r'[A-Z][0-9]{1,}',principal_hgvsp)
      if len(codon_match) == 1:
         principal_codon = codon_match[0]
   
   ## loop through all transcript-specific consequences ('csq_elements') for a given variant, check for the presence of
   ## 1. Exonic, protein-altering mutations (dictionary key = entrezgene + hgvsp) that overlap known cancer hotspots (https://github.com/sigven/cancerHotspots)
   ##    - assert whether a potentially hit reflects the principal hgvsp ('by_hgvsp_principal') or an alternative hgvsp ('by_hgvsp_nonprincipal')
   ## 2. Splice site mutations (dictionary key = entrezgene + hgvsc) that overlap known cancer hotspots (NOTE: only a subset have been curated)
   ##    - assert whether a potentially hit reflects the principal hgvsc ('by_hgvsc_principal') or an alternative hgvsc ('by_hgvsc_nonprincipal')

   for csq in transcript_csq_elements:
      (consequence, symbol, entrezgene, hgvsc, hgvsp, feature_type, feature, biotype) = csq.split(':')

      if not bool(re.search(r'^(missense|stop|start|inframe|splice_donor|splice_acceptor|frameshift)', consequence)) is True:
         continue

      hgvsp_short = threeToOneAA(hgvsp)
      hotspot_key_mutation = "."
      codon_match = []
      if entrezgene != "." and hgvsp != ".":
         hotspot_key_mutation = str(entrezgene) + '-' + str(hgvsp_short)
         codon_match = re.findall(r'p.[A-Z][0-9]{1,}',hgvsp_short)

      if entrezgene != "." and (consequence == 'splice_donor_variant' or consequence == 'splice_acceptor_variant'):
         hgvsc_key = re.sub(r'>(A|G|C|T)$','',hgvsc)
         hotspot_key_mutation = str(entrezgene) + '-' + str(hgvsc_key)

      if hotspot_key_mutation == ".":
         continue

      if hotspot_key_mutation in cancer_hotspots['mutation']:
         hotspot_info = cancer_hotspots['mutation'][hotspot_key_mutation]['MUTATION_HOTSPOT2']
         hotspot_info_ttype = cancer_hotspots['mutation'][hotspot_key_mutation]['MUTATION_HOTSPOT_CANCERTYPE']
         unique_hotspot_mutations['exonic|' + str(hotspot_info)] = hotspot_info_ttype

      if hotspot_key_mutation in cancer_hotspots['splice']:
         hotspot_info = cancer_hotspots['splice'][hotspot_key_mutation]['MUTATION_HOTSPOT2']
         hotspot_info_ttype = cancer_hotspots['splice'][hotspot_key_mutation]['MUTATION_HOTSPOT_CANCERTYPE']
         unique_hotspot_mutations['splice|' + str(hotspot_info)] = hotspot_info_ttype

               
      if len(codon_match) > 0:
         hotspot_key_codon = str(entrezgene) + '-' + str(codon_match[0])

         if hotspot_key_codon in cancer_hotspots['codon']:            
            unique_hotspot_codons[str('exonic|') + cancer_hotspots['codon'][hotspot_key_codon]['MUTATION_HOTSPOT2']] =  \
               cancer_hotspots['codon'][hotspot_key_codon]['MUTATION_HOTSPOT_CANCERTYPE']
         
   if len(unique_hotspot_mutations.keys()) > 0:
      if len(unique_hotspot_mutations.keys()) == 1:
         for gene_mutation_key in unique_hotspot_mutations.keys():
            rec.INFO['MUTATION_HOTSPOT'] = gene_mutation_key
            rec.INFO['MUTATION_HOTSPOT_CANCERTYPE'] = unique_hotspot_mutations[gene_mutation_key]

            if gene_mutation_key.split('|')[0] == 'exonic':
               rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsp_nonprincipal'
               hgvsp_candidate = 'p.' + str(gene_mutation_key.split('|')[3]) + str(gene_mutation_key.split('|')[4])
               if hgvsp_candidate == principal_hgvsp:
                  rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsp_principal'
            else:
               rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsc_nonprincipal'
               hgvsc_candidate = re.sub(r'>(A|G|C|T){1,}$', '' , str(gene_mutation_key.split('|')[4]))
               if hgvsc_candidate == principal_hgvsc:
                  rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsc_principal'
      else:
         ## multiple hotspot matches for alternative hgvsp keys
         ## - keep only match against principal HGVSp 
         for hotspot_info in unique_hotspot_mutations.keys():
            if hotspot_info.split('|')[0] == 'exonic':
               hgvsp_candidate = "p." + str(hotspot_info.split('|')[3]) + str(hotspot_info.split('|')[4]) 

               if hgvsp_candidate == principal_hgvsp:
                  rec.INFO['MUTATION_HOTSPOT'] = hotspot_info
                  rec.INFO['MUTATION_HOTSPOT_CANCERTYPE'] = unique_hotspot_mutations[hotspot_info]
                  rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsp_principal'
            else:
               hgvsc_candidate = re.sub(r'>(A|G|C|T){1,}$', '' , str(hotspot_info.split('|')[4]))

               if hgvsc_candidate == principal_hgvsc:
                  rec.INFO['MUTATION_HOTSPOT'] = hotspot_info
                  rec.INFO['MUTATION_HOTSPOT_CANCERTYPE'] = unique_hotspot_mutations[hotspot_info]
                  rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_hgvsc_principal'

   else:
      if len(unique_hotspot_codons.keys()) > 0:
         if len(unique_hotspot_codons.keys()) == 1:
            for gene_codon_key in unique_hotspot_codons.keys():

               if '|' in gene_codon_key:

                  codon = str(gene_codon_key.split('|')[3])

                  if codon == principal_codon:
                     rec.INFO['MUTATION_HOTSPOT'] = gene_codon_key
                     rec.INFO['MUTATION_HOTSPOT_CANCERTYPE'] = unique_hotspot_codons[gene_codon_key]
                     rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_codon_principal'
                  else:
                     rec.INFO['MUTATION_HOTSPOT'] = gene_codon_key
                     rec.INFO['MUTATION_HOTSPOT_CANCERTYPE'] = unique_hotspot_codons[gene_codon_key]
                     rec.INFO['MUTATION_HOTSPOT_MATCH'] = 'by_codon_nonprincipal'


   return
      

