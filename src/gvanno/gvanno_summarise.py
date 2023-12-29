#!/usr/bin/env python

import csv
import re
import argparse
import cyvcf2
import os

from lib.gvanno.annoutils import read_infotag_file, make_transcript_xref_map, read_genexref_namemap, write_pass_vcf, map_regulatory_variant_annotations
from lib.gvanno.vep import parse_vep_csq
from lib.gvanno.dbnsfp import vep_dbnsfp_meta_vcf, map_variant_effect_predictors
from lib.gvanno.oncogenicity import assign_oncogenicity_evidence
from lib.gvanno.mutation_hotspot import load_mutation_hotspots, match_csq_mutation_hotspot
from lib.gvanno.utils import error_message, check_subprocess, getlogger
from lib.gvanno.vep import parse_vep_csq

csv.field_size_limit(500 * 1024 * 1024)

def __main__():
    parser = argparse.ArgumentParser(description='Summarise VEP annotations (gene/variant) from gvanno pipeline (SNVs/InDels)')
    parser.add_argument('vcf_file_in', help='Bgzipped VCF file with VEP-annotated query variants (SNVs/InDels)')
    parser.add_argument('vcf_file_out', help='Bgzipped VCF file with extended VEP-annotated query variants (SNVs/InDels)')
    parser.add_argument('regulatory_annotation',default=0,type=int,help='Inclusion of VEP regulatory annotations (0/1)')
    parser.add_argument('oncogenicity_annotation',default=0,type=int,help='Include oncogenicity annotation (0/1)')
    parser.add_argument('vep_pick_order', default="mane_select,mane_plus_clinical,canonical,appris,biotype,ccds,rank,tsl,length", 
                        help=f"Comma-separated string of ordered transcript/variant properties for selection of primary variant consequence")
    parser.add_argument('gvanno_db_dir',help='gvanno data directory')
    parser.add_argument('--compress_output_vcf', action="store_true", default=False, help="Compress output VCF file")
    
    parser.add_argument("--debug", action="store_true", default=False, help="Print full commands to log, default: %(default)s")
    args = parser.parse_args()

    logger = getlogger('gvanno-gene-annotate')
    
    arg_dict = vars(args)
    
    extend_vcf_annotations(arg_dict, logger)

def extend_vcf_annotations(arg_dict, logger):
    """
    Function that reads VEP/vcfanno-annotated VCF and extends the VCF INFO column with tags from
    1. CSQ elements for the primary (i.e. "picked") gene transcript consequence from VEP, e.g. SYMBOL, Feature, Gene, Consequence etc.
    2. Gene annotations (GENE_TRANSCRIPT_XREF),
    3. Variant effect predictions - dbNSFP

    Moreover, it performs two important matching procedures, using 
    4. Information from VEP's CSQ information (HGVSp/HGVSc) to match known mutation hotspots in cancer
    
    Finally, it assesses somatic variant oncogenicity, using
    5. Gene annotations (tumor suppressor, oncogene) and variant annotations (loss-of-function, gnomAD variant frequencies, variant effect predictions).
       Variant oncogenicity levels are provided for all variants using a recommended five-level scheme ("Oncogenic", "Likely oncogenic", "VUS", "Likely Benign", "Benign")
       - Recommended scoring scheme for variant oncogenicity classification outlined by VICC/ClinGen consortia (Horak et al., Genet Med, 2022)

    List of VCF INFO tags appended by this procedure is defined by the 'infotags' files in the gvanno_db_dir
    """
    
    vcf_infotags = {}
    
    vcf_infotags['other'] = read_infotag_file(os.path.join(arg_dict['gvanno_db_dir'], 'vcf_infotags_gvanno.tsv'), scope = "gvanno")
    vcf_infotags['vep'] = read_infotag_file(os.path.join(arg_dict['gvanno_db_dir'], 'vcf_infotags_vep.tsv'), scope = "vep")
    vcf_infotags['other'].update(vcf_infotags['vep'])
    vcf_info_metadata = vcf_infotags['other']

    gene_transcript_xref_map = read_genexref_namemap(
        os.path.join(arg_dict['gvanno_db_dir'], 'gene','tsv','gene_transcript_xref', 'gene_transcript_xref_bedmap.tsv.gz'), logger)
    cancer_hotspots = load_mutation_hotspots(
        os.path.join(arg_dict['gvanno_db_dir'], 'misc','tsv','hotspot', 'hotspot.tsv.gz'), logger)

    out_vcf = re.sub(r'(\.gz)$','',arg_dict['vcf_file_out'])

    meta_vep_dbnsfp_info = vep_dbnsfp_meta_vcf(arg_dict['vcf_file_in'], vcf_info_metadata)
    dbnsfp_prediction_algorithms = meta_vep_dbnsfp_info['dbnsfp_prediction_algorithms']
    vep_csq_fields_map = meta_vep_dbnsfp_info['vep_csq_fieldmap']
    
    vcf = cyvcf2.VCF(arg_dict['vcf_file_in'])
    for tag in sorted(vcf_info_metadata):
        if arg_dict['regulatory_annotation'] == 0:
            if not tag.startswith('REGULATORY_'):
                vcf.add_info_to_header({'ID': tag, 'Description': str(vcf_info_metadata[tag]['description']),'Type':str(vcf_info_metadata[tag]['type']), 'Number': str(vcf_info_metadata[tag]['number'])})
        else:
            vcf.add_info_to_header({'ID': tag, 'Description': str(vcf_info_metadata[tag]['description']),'Type':str(vcf_info_metadata[tag]['type']), 'Number': str(vcf_info_metadata[tag]['number'])})

    w = cyvcf2.Writer(arg_dict['vcf_file_out'], vcf)
    current_chrom = None
    num_chromosome_records_processed = 0

    vcf_info_element_types = {}
    for e in vcf.header_iter():
        header_element = e.info()
        if 'ID' in header_element and 'HeaderType' in header_element and 'Type' in header_element:
            identifier = str(header_element['ID'])
            fieldtype = str(header_element['Type'])
            vcf_info_element_types[identifier] = fieldtype

    vars_no_csq = list()
    for rec in vcf:
        alt_allele = ','.join(rec.ALT)
        pos = rec.start + 1
        variant_id = f"g.{rec.CHROM}:{pos}{rec.REF}>{alt_allele}"
        if current_chrom is None:
            current_chrom = str(rec.CHROM)
            num_chromosome_records_processed = 0
        else:
            if str(rec.CHROM) != current_chrom:
                if not current_chrom is None:
                    logger.info(f"Completed summary of functional annotations for {num_chromosome_records_processed} variants on chr{current_chrom}")
                current_chrom = str(rec.CHROM)
                num_chromosome_records_processed = 0
        if rec.INFO.get('CSQ') is None:
            
            vars_no_csq.append(variant_id)
            continue

        num_chromosome_records_processed += 1
        transcript_xref_map = make_transcript_xref_map(rec, gene_transcript_xref_map, xref_tag = "GENE_TRANSCRIPT_XREF")

        vep_csq_record_results = {}
        vep_csq_record_results = \
            parse_vep_csq(rec, transcript_xref_map, vep_csq_fields_map, arg_dict['vep_pick_order'], 
                        logger, pick_only = False, csq_identifier = 'CSQ')
        
        if 'picked_gene_csq' in vep_csq_record_results and bool(arg_dict['regulatory_annotation']) is True:
            rec.INFO['REGULATORY_ANNOTATION'] = map_regulatory_variant_annotations(
                vep_csq_record_results['picked_gene_csq'])

        principal_csq_properties = {}
        principal_csq_properties['hgvsp'] = '.'
        principal_csq_properties['hgvsc'] = '.'
        principal_csq_properties['entrezgene'] = '.'
        principal_csq_properties['exon'] = '.'
        principal_csq_properties['codon'] = '.'
        principal_csq_properties['lof'] = '.'
        
        if 'picked_csq' in vep_csq_record_results:
            csq_record = vep_csq_record_results['picked_csq']
            for k in csq_record:
                if k in vcf_info_element_types:
                    if vcf_info_element_types[k] == "Flag" and csq_record[k] == "1":
                        rec.INFO[k] = True
                    else:
                        if not csq_record[k] is None:
                            rec.INFO[k] = csq_record[k]

                            if k == 'HGVSp_short':
                                principal_csq_properties['hgvsp'] = csq_record[k]
                                if re.match(r'^(p.[A-Z]{1}[0-9]{1,}[A-Za-z]{1,})', principal_csq_properties['hgvsp']):
                                    codon_match = re.findall(r'[A-Z][0-9]{1,}', principal_csq_properties['hgvsp'])
                                    if len(codon_match) == 1:
                                        principal_csq_properties['codon'] = 'p.' + codon_match[0]
                     
                            if k == 'HGVSc':
                                principal_csq_properties['hgvsc'] = csq_record[k].split(':')[1]
                            
                            if k == 'ENTREZGENE':
                                principal_csq_properties['entrezgene'] = csq_record[k]
                            
                            if k == 'LOSS_OF_FUNCTION':
                                principal_csq_properties['lof'] = csq_record[k]
                            
                            if k == 'EXON':
                                if "/" in csq_record[k]:
                                    principal_csq_properties['exon'] = csq_record[k].split('/')[0]
        
        if 'all_csq' in vep_csq_record_results:
            rec.INFO['VEP_ALL_CSQ'] = ','.join(vep_csq_record_results['all_csq'])
            match_csq_mutation_hotspot(vep_csq_record_results['all_csq'], cancer_hotspots, rec, principal_csq_properties)

        if not rec.INFO.get('DBNSFP') is None:
            map_variant_effect_predictors(rec, dbnsfp_prediction_algorithms)
        
        if arg_dict['oncogenicity_annotation'] == 1:
            assign_oncogenicity_evidence(rec, tumortype = "Any")

        if "GENE_TRANSCRIPT_XREF" in vcf_info_element_types:
            gene_xref_tag = rec.INFO.get('GENE_TRANSCRIPT_XREF')
            if not gene_xref_tag is None:
                del rec.INFO['GENE_TRANSCRIPT_XREF']                
        w.write_record(rec)
    if vars_no_csq:
        logger.warning(f"There were {len(vars_no_csq)} records with no CSQ tag from VEP (was --vep_no_intergenic flag set?). Skipping them and showing (up to) the first 100:")
        print('----')
        print(', '.join(vars_no_csq[:100]))
        print('----')
    w.close()
    if current_chrom is not None:
        logger.info(f"Completed summary of functional annotations for {num_chromosome_records_processed} variants on chr{current_chrom}")
    vcf.close()

    if os.path.exists(arg_dict['vcf_file_out']):
        if os.path.getsize(arg_dict['vcf_file_out']) > 0:
            if arg_dict['compress_output_vcf'] is True:
                check_subprocess(logger, f'bgzip -f {out_vcf}', debug=arg_dict['debug'])
                check_subprocess(logger, f'tabix -f -p vcf {out_vcf}.gz', debug=arg_dict['debug'])
                summarized_vcf = f'{arg_dict["vcf_file_out"]}.gz'
                write_pass_vcf(summarized_vcf, logger)
        else:
            error_message('No remaining PASS variants found in query VCF - exiting and skipping STEP 4', logger)
    else:
        error_message('No remaining PASS variants found in query VCF - exiting and skipping STEP 4', logger)

if __name__=="__main__":
    __main__()
