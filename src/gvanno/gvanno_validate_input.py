#!/usr/bin/env python

import csv
import re
import argparse
import os
import subprocess
import logging
import sys
import pandas as np
from cyvcf2 import VCF

from lib.gvanno.vcf import check_existing_vcf_info_tags
from lib.gvanno.annoutils import read_infotag_file
from lib.gvanno.utils import getlogger, random_id_generator, check_subprocess, remove_file


def __main__():
   
   parser = argparse.ArgumentParser(description='Verify input data for gvanno')
   parser.add_argument('gvanno_dir',help='Docker location of gvanno base directory with accompanying data directory, e.g. /data')
   parser.add_argument('input_vcf', help='VCF input file with query variants (SNVs/InDels)')
   parser.add_argument('validated_vcf', help="Validated VCF file with query variants (SNVs/InDels)")
   parser.add_argument('genome_assembly',help='grch37 or grch38')
   parser.add_argument('sample_id',help='Sample identifier')
   parser.add_argument('--output_dir', dest='output_dir', help='Output directory', default='/workdir/output')
   parser.add_argument('--debug', action='store_true', help="Print debug messages")
   args = parser.parse_args()
   
   ret = validate_gvanno_input(args.gvanno_dir, args.input_vcf, args.validated_vcf, args.sample_id, args.genome_assembly, args.output_dir, args.debug)
   if ret != 0:
      sys.exit(-1)



def simplify_vcf(input_vcf, validated_vcf, vcf, output_dir, sample_id, logger, debug):
    """
    input_vcf: path to input VCF
    validated_vcf: path to validated VCF
    vcf: parsed cyvcf2 object
    Function that performs the following on the validated input VCF:
    1. Strip of any genotype data
    2. If VCF has variants with multiple alternative alleles ("multiallelic", e.g. 'A,T'), 
       these are decomposed into variants with a single alternative allele
    3. Final VCF file is sorted and indexed (bgzip + tabix)
    """

    random_id = random_id_generator(15) 

    temp_files = {}
    temp_files['vcf_1'] = \
        os.path.join(output_dir, f'{sample_id}.gvanno_validate.bcftools.{random_id}_1.vcf')
    temp_files['vcf_2'] = \
        os.path.join(output_dir, f'{sample_id}.gvanno_validate.bcftools.{random_id}_2.vcf.gz')
    temp_files['vcf_3'] = \
        os.path.join(output_dir, f'{sample_id}.gvanno_validate.bftools.{random_id}_3.vcf.gz')
    bcftools_simplify_log = \
        os.path.join(output_dir, f'{sample_id}.gvanno_validate.bcftools.{random_id}.log')
    vt_decompose_log = \
        os.path.join(output_dir, f'{sample_id}.gvanno_validate.vt_decompose.{random_id}.log')

    multiallelic_list = list()
    for rec in vcf:
        POS = rec.start + 1
        alt = ",".join(str(n) for n in rec.ALT)
        if len(rec.ALT) > 1:
            variant_id = f"{rec.CHROM}:{POS}_{rec.REF}->{alt}"
            multiallelic_list.append(variant_id)

    logger.info('Extracting variants on autosomal/sex/mito chromosomes only (1-22,X,Y, M/MT) with bcftools')
    # bgzip + tabix required for sorting
    cmd_vcf1 = f'bcftools view {input_vcf} | bgzip -cf > {temp_files["vcf_2"]} && tabix -p vcf {temp_files["vcf_2"]} && ' + \
        f'bcftools sort --temp-dir {output_dir} -Oz {temp_files["vcf_2"]} > {temp_files["vcf_3"]} 2> {bcftools_simplify_log} && ' + \
        f'tabix -p vcf {temp_files["vcf_3"]}'
    # Keep only autosomal/sex/mito chrom (handle hg38 and hg19), remove FORMAT metadata lines, keep cols 1-8, sub chr prefix
    chrom_to_keep = [str(x) for x in [*range(1,23), 'X', 'Y', 'M', 'MT']]
    chrom_to_keep = ','.join([*['chr' + chrom for chrom in chrom_to_keep], *[chrom for chrom in chrom_to_keep]])
    cmd_vcf2 = f'bcftools view --regions {chrom_to_keep} {temp_files["vcf_3"]} | egrep -v \'^##FORMAT=\' ' + \
        f'| cut -f1-8 | sed \'s/^chr//\' > {temp_files["vcf_1"]}'

    check_subprocess(logger, cmd_vcf1, debug)
    check_subprocess(logger, cmd_vcf2, debug)

    if multiallelic_list:
        logger.warning(f"There were {len(multiallelic_list)} multiallelic sites detected. Showing (up to) the first 100:")
        print('----')
        print(', '.join(multiallelic_list[:100]))
        print('----')
        logger.info('Decomposing multi-allelic sites in input VCF file using \'vt decompose\'')
        command_decompose = f'vt decompose -s {temp_files["vcf_1"]} > {validated_vcf} 2> {vt_decompose_log}'
        check_subprocess(logger, command_decompose, debug)
    else:
        logger.info('All sites seem to be decomposed - skipping decomposition!')
        check_subprocess(logger, f'cp {temp_files["vcf_1"]} {validated_vcf}', debug)

    keep_uncompressed = False
    # need to keep uncompressed copy for vcf2maf.pl if selected
    bgzip_cmd = f"bgzip -cf {validated_vcf} > {validated_vcf}.gz" if keep_uncompressed else f"bgzip -f {validated_vcf}"
    check_subprocess(logger, bgzip_cmd, debug)
    check_subprocess(logger, f'tabix -p vcf {validated_vcf}.gz', debug)

    if os.path.exists(f'{validated_vcf}.gz') and os.path.getsize(f'{validated_vcf}.gz') > 0:
        vcf = VCF(f'{validated_vcf}.gz')
        i = 0
        for rec in vcf:
            i = i + 1
        if len(vcf.seqnames) == 0 or i == 0:
            logger.info('')
            logger.info("Input VCF contains NO valid variants after VCF cleaning - quitting workflow")
            logger.info('')
            exit(1)

    if not debug:
        remove_file(temp_files["vcf_1"])
        remove_file(temp_files["vcf_2"])
        remove_file(temp_files["vcf_3"])
        remove_file(temp_files["vcf_2"] + str('.tbi'))
        remove_file(temp_files["vcf_3"] + str('.tbi'))
        remove_file(bcftools_simplify_log)
        remove_file(vt_decompose_log)

def validate_gvanno_input(gvanno_directory, input_vcf, validated_vcf, sample_id, genome_assembly, output_dir, debug):
   """
   Function that reads the input file to gvanno (VCF file) and performs the following checks:
   1. Check that no INFO annotation tags in the query VCF coincides with those generated by gvanno
   2. Check that if VCF have variants with multiple alternative alleles (e.g. 'A,T') run vt decompose
   3. Any genotype data from VCF input file is stripped, and the resulting VCF file is sorted and indexed (bgzip + tabix) 
   """
   logger = getlogger('gvanno-validate-input')

   if not input_vcf == 'None':
      vcf_object = VCF(input_vcf)

      ## Check that VCF does not contain INFO tags that will be appended with gvanno annotation
      vcf_infotags = {}
      vcf_infotags['gvanno'] = read_infotag_file(
         os.path.join(gvanno_directory,'data',genome_assembly, 'vcf_infotags_gvanno.tsv'), scope = "gvanno")
      vcf_infotags['vep'] = read_infotag_file(
         os.path.join(gvanno_directory,'data',genome_assembly, 'vcf_infotags_vep.tsv'), scope = "vep")
      vcf_infotags['gvanno'].update(vcf_infotags['vep'])
      vcf_tags_gvanno = vcf_infotags['gvanno']
        
      tag_check = check_existing_vcf_info_tags(vcf_object, vcf_tags_gvanno, logger)
      if tag_check == -1:
         return -1     
      
      vcf = VCF(input_vcf)
      simplify_vcf(input_vcf, validated_vcf, vcf, output_dir, sample_id, logger, debug)
   
   return 0
   
if __name__=="__main__": __main__()

