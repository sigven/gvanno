#!/usr/bin/env python

import logging

from lib.gvanno.utils import error_message, warn_message, check_subprocess
from cyvcf2 import VCF
from typing import Union

def get_vcf_info_tags(vcf_fname):
    vcf = VCF(vcf_fname)
    info_tags = {}
    for e in vcf.header_iter():
        header_element = e.info()
        if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
            if header_element['HeaderType'] == 'INFO':
                info_tags[str(header_element['ID'])] = 1

    return info_tags


def print_vcf_header(vcf_fname, vcfheader_file, logger, chromline_only=False):
    if chromline_only == True:
        check_subprocess(
            logger, f'bgzip -dc {vcf_fname} | egrep \'^#\' | egrep \'^#CHROM\' >> {vcfheader_file}', debug=False)
    else:
        check_subprocess(
            logger, f'bgzip -dc {vcf_fname} | egrep \'^#\' | egrep -v \'^#CHROM\' > {vcfheader_file}', debug=False)

def detect_reserved_info_tag(tag, tag_name, logger):
    reserved_tags = ['AA', 'AC', 'AF', 'AN', 'BQ', 'CIGAR', 'DB', 'DP', 'END',
                     'H2', 'H3', 'MQ', 'MQ0', 'NS', 'SB', 'SOMATIC', 'VALIDATED', '1000G']
    if tag in reserved_tags:
        err_msg = f'Custom INFO tag ({tag_name}) needs another name - \'{tag}\' is a reserved field in the VCF specification (INFO)'
        return error_message(err_msg, logger)

    reserved_format_tags = ['GT', 'DP', 'FT', 'GL',
                            'GLE', 'GQ', 'PL', 'HQ', 'PS', 'PQ', 'EC', 'MQ']
    if tag in reserved_format_tags:
        err_msg = 'Custom INFO tag ({tag_name}) needs another name - \'{tag}\' is a reserved field in the VCF specification (FORMAT)'
        return error_message(err_msg, logger)

def check_retained_vcf_info_tags(vcf: VCF, retained_info_tags: str, logger: logging.Logger) -> int:

    """
    Function that compares the INFO tags in the query VCF and retained INFO tags set by the user as retained for output
    If any retained tag is not in query VCF, an error will be returned
    """

    tags = str(retained_info_tags).split(',')
    info_elements_query_vcf = []

    #vcf = VCF(input_vcf)
    logger.info('Checking if existing INFO tags of query VCF file matches retained INFO tags set by the user')
    ret = 1
    for e in vcf.header_iter():
        header_element = e.info()
        if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
            if header_element['HeaderType'] == 'INFO':
                info_elements_query_vcf.append(header_element['ID'])

    for t in tags:
        if not t in info_elements_query_vcf:
            err_msg = "Retained INFO tag '" + str(t) + "' not found among INFO tags in query VCF - make sure retained VCF INFO tags are set correctly"
            error_message(err_msg, logger)
        else:
            logger.info("Retained INFO tag '" + str(t) + "' detected among INFO tags in query VCF")
        
    return ret



def check_existing_vcf_info_tags(vcf: VCF, populated_infotags: dict, logger) -> int:

    """
    Function that compares the INFO tags in the query VCF and the INFO tags generated by gvanno
    If any coinciding tags, an error will be returned
    """

    logger.info('Checking if existing INFO tags of query VCF file coincide with gvanno INFO tags')
    ret = 1
    for e in vcf.header_iter():
        header_element = e.info()
        if 'ID' in header_element.keys() and 'HeaderType' in header_element.keys():
            if header_element['HeaderType'] == 'INFO':
                if header_element['ID'] in populated_infotags.keys():
                    err_msg = 'INFO tag ' + str(header_element['ID']) + ' in the query VCF coincides with a VCF ' + \
                        'annotation tag produced by gvanno - please remove or rename this tag in your query VCF'
                    return error_message(err_msg, logger)

    logger.info('No query VCF INFO tags coincide with populated INFO tags by gvanno')
    return ret
