#!/usr/bin/env python

GVANNO_VERSION = '1.7.0'
DB_VERSION = '20231222'

## MISCELLANEOUS
NCBI_BUILD_MAF = 'GRCh38'
MAX_VARIANTS_FOR_REPORT = 500_000
CODING_EXOME_SIZE_MB = 34.0
RECOMMENDED_N_MUT_SIGNATURE = 200

## GENCODE
GENCODE_VERSION = {'grch38': 44,'grch37': 19}

## vcfanno
VCFANNO_MAX_PROC = 15

## VEP settings/versions
VEP_VERSION = '110'
VEP_ASSEMBLY = {'grch38': 'GRCh38','grch37': 'GRCh37'}
VEP_MIN_FORKS = 1
VEP_MAX_FORKS = 8
VEP_MIN_BUFFER_SIZE = 50
VEP_MAX_BUFFER_SIZE = 30000
VEP_PICK_CRITERIA = ['mane_select','mane_plus_clinical','canonical','appris','tsl','biotype','ccds','rank','length']

## https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences
VEP_consequence_rank = {
    'transcript_ablation': 1,
    'splice_acceptor_variant': 2,
    'splice_donor_variant': 3,
    'stop_gained': 4,
    'frameshift_variant': 5,
    'stop_lost': 6,
    'start_lost' : 7,
    'transcript_amplification': 8,
    'feature_elongation': 9,
    'feature_truncation': 10,
    'inframe_insertion': 11,
    'inframe_deletion': 12,
    'missense_variant': 13,
    'protein_altering_variant': 14,
    'splice_donor_5th_base_variant': 15,
    'splice_region_variant': 16,
    'splice_donor_region_variant': 17,
    'splice_polypyrimidine_tract_variant': 18,
    'incomplete_terminal_codon_variant': 19,
    'start_retained_variant': 20,
    'stop_retained_variant': 21,
    'synonymous_variant': 22,
    'coding_sequence_variant': 23,
    'mature_miRNA_variant': 24,
    '5_prime_UTR_variant': 25,
    '3_prime_UTR_variant': 26,
    'non_coding_transcript_exon_variant': 27,
    'intron_variant': 28,
    'NMD_transcript_variant': 29,
    'non_coding_transcript_variant': 30,
    'coding_transcript_variant': 31,
    'upstream_gene_variant': 32,
    'downstream_gene_variant': 33,
    'TFBS_ablation': 34,
    'TFBS_amplification': 35,
    'TF_binding_site_variant': 36,
    'regulatory_region_ablation': 37,
    'regulatory_region_amplification': 38,
    'regulatory_region_variant': 39,
    'intergenic_variant': 40,
    'sequence_variant': 41
}

CSQ_MISSENSE_PATTERN = r"^missense_variant"
CSQ_CODING_PATTERN = r"(stop_(lost|gained)|start_lost|frameshift_|missense_|splice_(donor|acceptor)|protein_altering|inframe_)"
CSQ_CODING_SILENT_PATTERN = r"(stop_(lost|gained)|start_lost|frameshift_|missense_|splice_(donor|acceptor)|protein_altering|inframe_|synonymous|(start|stop)_retained)"
CSQ_NULL_PATTERN = r"(stop_gained|frameshift_)"
CSQ_SPLICE_REGION_PATTERN = r"(splice_|intron_variant)"
CSQ_SPLICE_DONOR_PATTERN = r"(splice_region_variant|splice_donor_variant|splice_donor_region_variant|splice_donor_5th_base_variant)"

