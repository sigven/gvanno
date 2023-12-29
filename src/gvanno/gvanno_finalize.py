#!/usr/bin/env python

import os,re
import argparse

from lib.gvanno.variant import append_annotations, clean_annotations
from lib.gvanno.utils import getlogger


def __main__():
    parser = argparse.ArgumentParser(description='Finalize gvanno-annotated TSV file - append final annotations')
    parser.add_argument('gvanno_db_dir', help='Assembly-specific directory with gvanno data bundle files')
    parser.add_argument('tsv_file_in', help='vcf2tsv-converted VCF file with gvanno-annotated variants (SNVs/InDels)')
    parser.add_argument('tsv_file_out', help='TSV file with cleaned gvanno-annotated variants (SNVs/InDels)')
    parser.add_argument('genome_assembly', help='Genome assembly')
    parser.add_argument('sample_id', help='Sample identifier')
    parser.add_argument("--debug", action="store_true", default=False, help="Print full commands to log, default: %(default)s")
    args = parser.parse_args()

    logger = getlogger('gvanno-vep')
    
    arg_dict = vars(args)
   
    variant_set = \
        append_annotations(
            arg_dict['tsv_file_in'], gvanno_db_dir = arg_dict['gvanno_db_dir'], logger = logger)
    variant_set = clean_annotations(variant_set, arg_dict['sample_id'], arg_dict['genome_assembly'], logger = logger)        
    variant_set.fillna('.').to_csv(arg_dict['tsv_file_out'], sep="\t", compression="gzip", index=False)
    

if __name__=="__main__": __main__()

