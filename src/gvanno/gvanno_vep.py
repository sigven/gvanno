#!/usr/bin/env python

import os,re
import argparse

from lib.gvanno.utils import get_loftee_dir, getlogger, check_subprocess
from lib.gvanno import gvanno_vars


def __main__():
    parser = argparse.ArgumentParser(description='Run VEP on a VCF file (SNVs/InDels)')
    parser.add_argument('vep_cache_dir', help='Directory with VEP cache files)')
    parser.add_argument('vcf_file_in', help='Bgzipped VCF file with validated query variants (SNVs/InDels)')
    parser.add_argument('vcf_file_out', help='Bgzipped VCF file with VEP-annotated query variants (SNVs/InDels)')
    parser.add_argument('genome_assembly', help='Genome assembly')
    parser.add_argument('vep_pick_order', default="mane_select,mane_plus_clinical,canonical,appris,biotype,ccds,rank,tsl,length", 
                        help=f"Comma-separated string of ordered transcript/variant properties for selection of primary variant consequence")
    parser.add_argument('--vep_regulatory',action = "store_true",help='Inclusion of VEP regulatory annotations')
    
    parser.add_argument('--vep_buffer_size',default=1000,type=int,help='Buffer size for VEP')
    parser.add_argument('--vep_gencode_basic',action="store_true",help='Only consier basic GENCODE transcripts)')
    parser.add_argument('--vep_n_forks',default=4,type=int,help='Number of forks for VEP processing')
    parser.add_argument('--vep_lof_prediction',action="store_true",help='Perform LoF prediction with the LOFTEE plugin in VEP')
    parser.add_argument('--vep_coding_only', action="store_true", help="Only consider coding variants")
    parser.add_argument('--vep_no_intergenic', action="store_true", help="Skip intergenic variants")
    parser.add_argument("--debug", action="store_true", default=False, help="Print full commands to log, default: %(default)s")
    args = parser.parse_args()

    logger = getlogger('gvanno-vep')
    
    arg_dict = vars(args)
    conf_options = {}
    conf_options['genome_assembly'] = arg_dict['genome_assembly']
    conf_options['conf'] = {}
    conf_options['conf']['vep'] = {}
    conf_options['conf']['vep']['vep_n_forks'] = arg_dict['vep_n_forks']
    conf_options['conf']['vep']['vep_pick_order'] = arg_dict['vep_pick_order']
    conf_options['conf']['vep']['vep_buffer_size'] = arg_dict['vep_buffer_size']
    conf_options['conf']['vep']['vep_gencode_basic'] = arg_dict['vep_gencode_basic']
    conf_options['conf']['vep']['vep_regulatory'] = arg_dict['vep_regulatory']
    conf_options['conf']['vep']['vep_lof_prediction'] = arg_dict['vep_lof_prediction']
    conf_options['conf']['vep']['vep_coding_only'] = arg_dict['vep_coding_only']
    conf_options['conf']['vep']['vep_no_intergenic'] = arg_dict['vep_no_intergenic']
    
    run_vep(arg_dict['vep_cache_dir'], conf_options, arg_dict['vcf_file_in'], arg_dict['vcf_file_out'], logger, arg_dict['debug'])

def run_vep(vep_cache_dir, conf_options, input_vcf, output_vcf, logger, debug = False):
    
    output_vcf_gz = f'{output_vcf}.gz'
    genome_assembly = conf_options['genome_assembly']
    
    fasta_assembly = os.path.join(
        vep_cache_dir, 'homo_sapiens', 
        f'{gvanno_vars.VEP_VERSION}_{gvanno_vars.VEP_ASSEMBLY[genome_assembly]}', 
        f'Homo_sapiens.{gvanno_vars.VEP_ASSEMBLY[genome_assembly]}.dna.primary_assembly.fa.gz')
    ancestor_assembly = os.path.join(
        vep_cache_dir, 'homo_sapiens', 
        f'{gvanno_vars.VEP_VERSION}_{gvanno_vars.VEP_ASSEMBLY[genome_assembly]}', 
        f'human_ancestor.fa.gz')

    plugins_in_use = "NearestExonJB, LoF"

    # List all VEP flags used when calling VEP
    vep_flags = (
        f'--hgvs --af_gnomad --variant_class --domains --symbol --protein --ccds --mane '
        f'--uniprot --appris --biotype --tsl --canonical --format vcf --cache --numbers '
        f'--total_length --allele_number --failed 1 --no_stats --no_escape --xref_refseq --vcf '
        f'--check_ref --dont_skip --flag_pick_allele_gene --plugin NearestExonJB,max_range=50000 '
        f'--force_overwrite --species homo_sapiens --offline')
    vep_options = (
        f'--dir {vep_cache_dir} --assembly {gvanno_vars.VEP_ASSEMBLY[genome_assembly]} --cache_version {gvanno_vars.VEP_VERSION} '
        f'--fasta {fasta_assembly} '
        f'--pick_order {conf_options["conf"]["vep"]["vep_pick_order"]} '
        f'--buffer_size {conf_options["conf"]["vep"]["vep_buffer_size"]} '
        f'--fork {conf_options["conf"]["vep"]["vep_n_forks"]} '
        f'{vep_flags} '
        f'{"--verbose" if debug else "--quiet"}')
    
    gencode_set_in_use = "GENCODE - all transcripts"
    if conf_options['conf']['vep']['vep_no_intergenic'] is True:
        vep_options += ' --no_intergenic'
    if conf_options['conf']['vep']['vep_regulatory'] is True:
        vep_options += ' --regulatory'
    if conf_options['conf']['vep']['vep_coding_only'] is True:
        vep_options += ' --coding_only'
    if conf_options['conf']['vep']['vep_gencode_basic'] is True:
        vep_options += ' --gencode_basic'
        gencode_set_in_use = "GENCODE - basic transcript set (--gencode_basic)"

    ## LOFTEE plugin - variant loss-of-function annotation        
    loftee_dir = get_loftee_dir()
    assert os.path.isdir(loftee_dir), f'LoF VEP plugin is not found in {loftee_dir}. Please make sure you installed pcgr conda package and have corresponding conda environment active.'
    vep_options += f" --plugin LoF,loftee_path:{loftee_dir},human_ancestor_fa:{ancestor_assembly},use_gerp_end_trunc:0 --dir_plugins {loftee_dir}"

    logger.info(f'VEP configuration: Version: {gvanno_vars.VEP_VERSION}' + \
                f', GENCODE release {gvanno_vars.GENCODE_VERSION[genome_assembly]}, genome assembly {conf_options["genome_assembly"]}')
    logger.info(f'VEP configuration - one primary consequence block pr. alternative allele (--flag_pick_allele)')
    logger.info(f'VEP configuration - transcript pick order: {conf_options["conf"]["vep"]["vep_pick_order"]}')
    logger.info(f'VEP configuration - transcript pick order: See more at https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options')
    logger.info(f'VEP configuration - GENCODE set: {gencode_set_in_use}')
    logger.info(f'VEP configuration - skip intergenic variants: {"ON" if conf_options["conf"]["vep"]["vep_no_intergenic"] == 1 else "OFF"}')
    logger.info(f'VEP configuration - regulatory variant annotation: {"ON" if conf_options["conf"]["vep"]["vep_regulatory"] == 1 else "OFF"}')
    logger.info(f'VEP configuration - loss-of-function prediction: {"ON" if conf_options["conf"]["vep"]["vep_lof_prediction"] == 1 else "OFF"}')

    logger.info((
        f'VEP configuration - buffer size/number of forks: '
        f'{conf_options["conf"]["vep"]["vep_buffer_size"]}/{conf_options["conf"]["vep"]["vep_n_forks"]}'))
    logger.info(f'VEP - plugins in use: {plugins_in_use}')
    
    # Compose full VEP command
    vep_main_command = f'vep --input_file {input_vcf} --output_file {output_vcf} {vep_options}'
    vep_bgzip_command = f'bgzip -f -c {output_vcf} > {output_vcf_gz}'
    vep_tabix_command = f'tabix -f -p vcf {output_vcf_gz}'
    if debug:
        print(vep_main_command)
    
    check_subprocess(logger, vep_main_command, debug)
    check_subprocess(logger, vep_bgzip_command, debug)
    check_subprocess(logger, vep_tabix_command, debug)
    logger.info('Finished gvanno-vep')
    
    return 0

if __name__=="__main__": __main__()

