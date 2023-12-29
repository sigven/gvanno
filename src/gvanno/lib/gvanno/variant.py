#!/usr/bin/env python

import os
import re
import pandas as pd
import numpy as np
import warnings

warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


def set_genotype(variant_set: pd.DataFrame, logger) -> pd.DataFrame:
    """
    Set verbose genotype (homozygous, heterozygous) for each variant
    """
    variant_set['GENOTYPE'] = '.'    
    if {'GT'}.issubset(variant_set.columns):
        logger.info("Assignment of genotype (homozygous, heterozygous) for each variant based on 'GT' tag")
        heterozygous_states = []
        ref_allele_index = 0
        while ref_allele_index < 20:
            alt_allele_index = ref_allele_index + 1
            while alt_allele_index <= 20:
                phased_gt_1 = str(ref_allele_index) + "|" + str(alt_allele_index)
                phased_gt_2 = str(alt_allele_index) + "|" + str(ref_allele_index)
                unphased_gt_1 = str(ref_allele_index) + "/" + str(alt_allele_index)
                unphased_gt_2 = str(alt_allele_index) + "/" + str(ref_allele_index)
                gts = [phased_gt_1, phased_gt_2, unphased_gt_1, unphased_gt_2]
                heterozygous_states.extend(gts)
                
                alt_allele_index = alt_allele_index + 1
            ref_allele_index = ref_allele_index + 1
            
        homozygous_states = []    
        hom_allele_index = 1
        while hom_allele_index <= 20:
            phased_gt = str(hom_allele_index) + "|" + str(hom_allele_index)
            unphased_gt = str(hom_allele_index) + "/" + str(hom_allele_index)
            homozygous_states.extend([phased_gt, unphased_gt])
            hom_allele_index = hom_allele_index + 1
        
        variant_set.loc[variant_set['GT'].isin(homozygous_states), "GENOTYPE"] = "homozygous"
        variant_set.loc[variant_set['GT'].isin(heterozygous_states),"GENOTYPE"] = "heterozygous"
    else:
        variant_set['GENOTYPE'] = "undefined"
    
    variant_set = variant_set.astype({'GENOTYPE':'string'})  
    return(variant_set)

def append_annotations(vcf2tsv_gz_fname: str, gvanno_db_dir: str, logger):
    
    clinvar_tsv_fname = os.path.join(gvanno_db_dir, 'variant','tsv','clinvar', 'clinvar.tsv.gz')
    protein_domain_tsv_fname = os.path.join(gvanno_db_dir, 'misc','tsv','protein_domain', 'protein_domain.tsv.gz') 
    gene_xref_tsv_fname = os.path.join(gvanno_db_dir, 'gene','tsv','gene_transcript_xref', 'gene_transcript_xref.tsv.gz')
    vcf2tsv_df = None
    clinvar_data_df = None
    
    num_recs_with_clinvar_hits = 0
    
    if os.path.exists(vcf2tsv_gz_fname):
        vcf2tsv_df = pd.read_csv(
            vcf2tsv_gz_fname, skiprows=[0], sep="\t", na_values=".",
            low_memory = False)
        if {'CHROM','POS','REF','ALT','CLINVAR_MSID','PFAM_DOMAIN','ENTREZGENE'}.issubset(vcf2tsv_df.columns):
            for elem in ['CHROM','POS','REF','ALT','CLINVAR_MSID','PFAM_DOMAIN','ENTREZGENE']:
                vcf2tsv_df = vcf2tsv_df.astype({elem:'string'})            
            vcf2tsv_df['CLINVAR_MSID'] = vcf2tsv_df['CLINVAR_MSID'].str.replace("\\.[0-9]{1,}$", "", regex = True)
            vcf2tsv_df['PFAM_DOMAIN'] = vcf2tsv_df['PFAM_DOMAIN'].str.replace("\\.[0-9]{1,}$", "", regex = True)
            vcf2tsv_df['ENTREZGENE'] = vcf2tsv_df['ENTREZGENE'].str.replace("\\.[0-9]{1,}$", "", regex = True)
            vcf2tsv_df["VAR_ID"] = vcf2tsv_df["CHROM"].str.cat(
                vcf2tsv_df["POS"], sep = "_").str.cat(
                    vcf2tsv_df["REF"], sep = "_").str.cat(
                        vcf2tsv_df["ALT"], sep = "_")

            if {'CLINVAR_TRAITS_ALL'}.issubset(vcf2tsv_df.columns):
                vcf2tsv_df.drop('CLINVAR_TRAITS_ALL', inplace=True, axis=1)
        
            ## check number of variants with ClinVar ID's
            num_recs_with_clinvar_hits = vcf2tsv_df["CLINVAR_MSID"].notna().sum()
            ## check number of variants with PFAM ID's
            num_recs_with_pfam_hits = vcf2tsv_df["PFAM_DOMAIN"].notna().sum()
            ## check number of variants with Ensembl gene ID's
            num_recs_with_entrez_hits = vcf2tsv_df["ENTREZGENE"].notna().sum()
    
            #print(str(num_recs_with_entrez_hits))
            ## merge variant set with ClinVar trait and variant origin annotations
            if num_recs_with_clinvar_hits > 0:
                if os.path.exists(clinvar_tsv_fname):
                    clinvar_data_df = pd.read_csv(
                        clinvar_tsv_fname, sep="\t", 
                        usecols=["variation_id","origin_simple","VAR_ID","trait"],
                        low_memory = False)
                    clinvar_data_df['CLINVAR_TRAITS_ALL'] = clinvar_data_df['origin_simple'].str.capitalize().str.cat(
                        clinvar_data_df['trait'], sep = " - ")
                    clinvar_data_df['CLINVAR_MSID'] = clinvar_data_df['variation_id']
                    clinvar_data_df = clinvar_data_df.astype({'CLINVAR_MSID':'string'})
                    clinvar_data_df['CLINVAR_MSID'] = clinvar_data_df['CLINVAR_MSID'].str.replace("\\.[0-9]{1,}$", "", regex = True)
                    clinvar_data_df = clinvar_data_df[['VAR_ID','CLINVAR_MSID','CLINVAR_TRAITS_ALL']]
                    
                    vcf2tsv_df = vcf2tsv_df.merge(
                        clinvar_data_df, left_on=["VAR_ID", "CLINVAR_MSID"], right_on=["VAR_ID", "CLINVAR_MSID"], how="left")
                else:
                    logger.error(f"Could not find {clinvar_tsv_fname} needed for ClinVar variant annotation - exiting")
            else:
                vcf2tsv_df['CLINVAR_TRAITS_ALL'] = '.'
                
            
            ## merge variant set with PFAM domain annotations
            if num_recs_with_pfam_hits > 0:
                
                vcf2tsv_df.drop('PFAM_DOMAIN_NAME', inplace=True, axis=1)
                
                if os.path.exists(protein_domain_tsv_fname):
                    prot_domains_data_df = pd.read_csv(
                        protein_domain_tsv_fname, sep="\t", usecols=["pfam_id","pfam_name"]).drop_duplicates()
                    prot_domains_data_df.rename(columns = {'pfam_id':'PFAM_DOMAIN', 'pfam_name':'PFAM_DOMAIN_NAME'}, inplace = True)                                        
                    vcf2tsv_df = vcf2tsv_df.merge(prot_domains_data_df, left_on=["PFAM_DOMAIN"], right_on=["PFAM_DOMAIN"], how="left")
                else:
                    logger.error(f"Could not find {protein_domain_tsv_fname} needed for PFAM domain annotation - exiting")
            else:
                vcf2tsv_df['PFAM_DOMAIN_NAME'] = '.'
            
            if num_recs_with_entrez_hits > 0:
                
                if {'GENENAME'}.issubset(vcf2tsv_df.columns):
                    vcf2tsv_df.drop('GENENAME', inplace=True, axis=1)
                
                if os.path.exists(gene_xref_tsv_fname):
                    gene_xref_df = pd.read_csv(
                        gene_xref_tsv_fname, sep="\t", na_values=".", 
                        usecols=["entrezgene","name"])
                    gene_xref_df = gene_xref_df[gene_xref_df['entrezgene'].notnull()].drop_duplicates()
                    gene_xref_df = gene_xref_df[gene_xref_df['entrezgene'].notna()].drop_duplicates()
                    gene_xref_df["entrezgene"] = gene_xref_df["entrezgene"].astype(float).astype(int).astype(str)
                    vcf2tsv_df["ENTREZGENE"] = vcf2tsv_df["ENTREZGENE"].astype(str)
                    vcf2tsv_df.loc[vcf2tsv_df["ENTREZGENE"].isna(), "ENTREZGENE"] = "-1"
                    gene_xref_df.rename(columns = {'entrezgene':'ENTREZGENE', 'name':'GENENAME'}, inplace = True)                
                    vcf2tsv_df = vcf2tsv_df.merge(gene_xref_df, left_on=["ENTREZGENE"], right_on=["ENTREZGENE"], how="left")
                    vcf2tsv_df["ENTREZGENE"] = vcf2tsv_df['ENTREZGENE'].str.replace("\\.[0-9]{1,}$", "", regex = True)
                else:
                    logger.error(f"Could not find {gene_xref_tsv_fname} needed for gene name annotation - exiting")
            else:
                vcf2tsv_df['GENENAME'] = '.'
    
    return(vcf2tsv_df)

def clean_annotations(variant_set: pd.DataFrame, sample_id, genome_assembly, logger) -> pd.DataFrame:
    
    if {'Consequence','EFFECT_PREDICTIONS','CLINVAR_CONFLICTED'}.issubset(variant_set.columns):
        variant_set.rename(columns = {'Consequence':'CONSEQUENCE'}, inplace = True)
        variant_set['EFFECT_PREDICTIONS'] = variant_set['EFFECT_PREDICTIONS'].str.replace("\\.&|\\.$", "NA&", regex = True)
        variant_set['EFFECT_PREDICTIONS'] = variant_set['EFFECT_PREDICTIONS'].str.replace("&$", "", regex = True)
        variant_set['EFFECT_PREDICTIONS'] = variant_set['EFFECT_PREDICTIONS'].str.replace("&", ", ", regex = True)
        variant_set['clinvar_conflicted_bool'] = True
        variant_set.loc[variant_set['CLINVAR_CONFLICTED'] == 1, "clinvar_conflicted_bool"] = True
        variant_set.loc[variant_set['CLINVAR_CONFLICTED'] != 1, "clinvar_conflicted_bool"] = False
        variant_set.drop('CLINVAR_CONFLICTED', inplace=True, axis=1)        
        variant_set.rename(columns = {'clinvar_conflicted_bool':'CLINVAR_CONFLICTED'}, inplace = True)
        
    
    if not {'VCF_SAMPLE_ID'}.issubset(variant_set.columns):
        variant_set['VCF_SAMPLE_ID'] = str(sample_id)
    variant_set['SAMPLE_ID'] = str(sample_id)
    variant_set['GENOME_VERSION'] = str(genome_assembly)
    if {'CHROM','POS','REF','ALT',}.issubset(variant_set.columns):      
        variant_set['GENOMIC_CHANGE'] = variant_set['CHROM'].astype(str) + ":g." + variant_set['POS'].astype(str) + \
            variant_set['REF'].astype(str) + ">" + variant_set['ALT'].astype(str)
   
    ## Make sure that specific tags are formatted as integers (not float) during to_csv export
    if {'AMINO_ACID_END','AMINO_ACID_START'}.issubset(variant_set.columns):
        variant_set.loc[variant_set['AMINO_ACID_START'].isna(),"AMINO_ACID_START"] = -1        
        variant_set.loc[variant_set['AMINO_ACID_END'].isna(),"AMINO_ACID_END"] = -1
        variant_set['AMINO_ACID_END'] = variant_set['AMINO_ACID_END'].astype(float).astype(int)
        variant_set['AMINO_ACID_START'] = variant_set['AMINO_ACID_START'].astype(float).astype(int)
    
    for elem in ['NUM_SUBMITTERS','ALLELE_ID','ENTREZGENE','REVIEW_STATUS_STARS']:
        vcf_info_tag = 'CLINVAR_' + str(elem)
        if vcf_info_tag in variant_set.columns:
            variant_set.loc[variant_set[vcf_info_tag].notna(), vcf_info_tag] = \
                variant_set.loc[variant_set[vcf_info_tag].notna(), vcf_info_tag].astype(str).astype(float).astype(int)
    
    for vcf_info_tag in ['ONCOGENE_RANK','TSG_RANK','TCGA_PANCANCER_COUNT','CGC_TIER','DISTANCE',
                         'EXON_AFFECTED','INTRON_POSITION','EXON_POSITION']:
        if vcf_info_tag in variant_set.columns:                        
            variant_set.loc[variant_set[vcf_info_tag].notna(), vcf_info_tag] = \
                variant_set.loc[variant_set[vcf_info_tag].notna(), vcf_info_tag].astype(str).astype(float).astype(int)
    
    variant_set = set_genotype(variant_set, logger)
    
    return variant_set