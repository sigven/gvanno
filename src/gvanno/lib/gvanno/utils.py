#!/usr/bin/env python
import sys
import subprocess
import shutil
import logging
import os
import time
import errno
import platform
import string
import random


def getlogger(logger_name):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to debug
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    # add ch to logger
    logger.addHandler(ch)
    # create formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s", "20%y-%m-%d %H:%M:%S")
    # add formatter to ch
    ch.setFormatter(formatter)
    return logger

def error_message(message, logger):
    logger.error("")
    logger.error(message)
    logger.error("")
    sys.exit(1)

def warn_message(message, logger):
    logger.warning("")
    logger.warning(message)
    logger.warning("")
 
def random_id_generator(size = 10, chars = string.ascii_lowercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))   

#def random_string(length):
#    letters = string.ascii_lowercase
#    result_str = ''.join(random.choice(letters) for i in range(length))
#    return result_str

def check_subprocess(logger, command, debug):
    if debug:
        logger.info(command)
    try:
        output = subprocess.check_output(
            str(command), stderr=subprocess.STDOUT, shell=True)
        if len(output) > 0:
            print(str(output.decode()).rstrip())
    except subprocess.CalledProcessError as e:
        print(e.output.decode())
        exit(0)

def get_loftee_dir():
    #pcgr_conda_env = conda_prefix_basename()
    return "/opt/vep/src/ensembl-vep/modules"

def perl_cmd():
    """Retrieve abs path to locally installed conda Perl or first in PATH,
    e.g. conda/env/pcgr/bin/perl
    """
    perl = shutil.which(os.path.join(get_pcgr_bin(), "perl"))
    if perl:
        return perl
    else:
        return shutil.which("perl")

def get_perl_exports():
    """Environmental exports to use conda installed perl.
    """
    perl_path = os.path.normpath(perl_cmd())      # /conda/env/pcgr/bin/perl
    perl_path_parent = os.path.dirname(perl_path) # /conda/env/pcgr/bin
    out = f"unset PERL5LIB && export PATH={perl_path_parent}:\"$PATH\""
    return out

def is_integer(n):
    try:
        float(n)
    except ValueError:
        return False
    else:
        return float(n).is_integer()

# https://stackoverflow.com/a/10840586/2169986
def remove_file(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred


# from https://github.com/bcbio/bcbio-nextgen/blob/master/bcbio/utils.py
def safe_makedir(dname):
    """Make a directory if it does not exist, handling concurrent race conditions.
    """
    if not dname:
        return dname
    num_tries = 0
    max_tries = 5
    while not os.path.exists(dname):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dname)
        except OSError:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dname

def sort_bed(unsorted_bed_fname: str, sorted_bed_fname: str, debug = False, logger = None):
    
    ## Sort regions in target BED
    if os.path.exists(unsorted_bed_fname) and os.stat(unsorted_bed_fname).st_size != 0:
        cmd_sort_custom_bed1 = 'egrep \'^[0-9]\' ' + str(unsorted_bed_fname) + \
            ' | sort -k1,1n -k2,2n -k3,3n > ' + str(sorted_bed_fname)
        cmd_sort_custom_bed2 = 'egrep -v \'^[0-9]\' ' + str(unsorted_bed_fname) + \
            ' | egrep \'^[XYM]\' | sort -k1,1 -k2,2n -k3,3n >> ' + str(sorted_bed_fname)

        check_subprocess(logger, cmd_sort_custom_bed1, debug)
        check_subprocess(logger, cmd_sort_custom_bed2, debug)
        if not debug:
            remove(str(unsorted_bed_fname))
    else:
        err_msg = 'File ' + str(unsorted_bed_fname) + ' does not exist or is empty'
        error_message(err_msg, logger)


def check_file_exists(fname: str, logger = None) -> bool:
    err = 0
    if not os.path.isfile(fname):
        err = 1
    else:
        if os.stat(fname).st_size == 0:
            err = 1
    if err == 1:
        err_msg = f"File {fname} does not exist or has zero size"
        error_message(err_msg, logger)
    
    return(True)

def check_tabix_file(fname: str, logger = None) -> bool:
    tabix_file = os.path.abspath(f'{fname}.tbi')
    if not os.path.exists(tabix_file):
        err_msg = f"Tabix file (i.e. '.gz.tbi') is not present for the bgzipped VCF input file ({fname}" + \
            "). Please make sure your input VCF is properly compressed and indexed (bgzip + tabix)"
        error_message(err_msg, logger)
    else:
        ## check file size is more than zero
        check_file_exists(tabix_file)
    return(True)
