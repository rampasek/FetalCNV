#!/usr/bin/python2

import argparse
import random
import copy
import samParse as sp

def main():
    #parse ARGs
    parser = argparse.ArgumentParser(description='Prepare SNP data for FCNV. Read filenames for M, P, and F .vcf files and plasma, M, and P .sam files.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to .vcf files with M, P, F SNPs and reads in SAM format for plasma, M, and P samples.')
    args = parser.parse_args()
    
    if len(args.filenames) != 6: die("Unexpected number of arguments passed! Expecting 6 filenames.")
    
    #treat these as CONSTANTS!
    M = 0; P = 1; F = 2; PLASMA = 3; MR = 4; PR = 5;
    ALL = [M, P, F, PLASMA, MR, PR]
    
    #list of input files
    in_files = [open(args.filenames[i], "r" ) for i in ALL]
    
    #list of output files
    out_files = [None for i in [M, P, F, PLASMA]]
    out_files[M] = open("M_alleles.txt", "w")
    out_files[P] = open("P_alleles.txt", "w")
    out_files[F] = open("F_alleles.txt", "w")
    out_files[PLASMA] = open("plasma_samples.txt", "w")
    out_pos_file = open("positions.txt", "w")
    
    #allele counts in plasma samples for particular positions
    data = dict()
    loci = dict()
    
    #read SNPs from M, P, F vcf files
    snps = [[] for i in [M, P, F]]
    for i in [M, P, F]:
        #skip the header
        line = in_files[i].readline()
        while len(line) > 0 and line[0] == '#': line = in_files[i].readline()
        #split
        snps[i] = line.split('\t')
    
    #output genotypes for all positions in UNION of M and P SNP positions
    while len(snps[M])>2 or len(snps[P])>2: #while there is a SNP positions in M or P
        #get the position of SNP that occure first
        for i in [M, P, F]:
            if snps[i][0] == '': #if an input files is already at EOF
                snps[i][0] = 'chrZZ'
                snps[i].append(1e15)
            else: #convert to int
                snps[i][1] = int(snps[i][1])
        #chromosome
        min_chr = min(snps[M][0], snps[P][0])
        #position
        min_pos = 1e15
        for i in [M, P]:
            if min_chr == snps[i][0] and snps[i][1] < min_pos:
                min_pos = snps[i][1]
        
        #get alleles
        alleles = [['.', '.'] for x in [M, P, F]]
        for i in [M, P, F]:
            #if there is a SNP in the data at this position, use it
            if min_chr == snps[i][0] and min_pos == snps[i][1]:
                alleles[i] = [snps[i][3], snps[i][4]]

            
        #if there is no info on some allele, impute the reference allele
        ref_allele = alleles[M][0]
        if ref_allele == '.': ref_allele = alleles[P][0]
        for i in [M, P, F]:
            for a in [0, 1]:
                if alleles[i][a] == '.': alleles[i][a] = ref_allele

        #check for homozygous alternative sites in Fetal VCF
        for i in [F]:
            #if there is a SNP in the data at this position
            if min_chr == snps[i][0] and min_pos == snps[i][1]:
                info = snps[i]
                if len(info) <= 2: continue
                #parse out genotype config info
                gt = info[9].split(':')[0]
                #if homozygous alternative
                if gt[0] == '1':
                    alleles[i][0] = snps[i][4]
                if gt[2] == '1':
                    alleles[i][1] = snps[i][4]
                    
        #organize the haplotypes in M, P (phased VCF files)
        for i in [M, P]:
            #if there is a SNP in the data at this position
            if min_chr == snps[i][0] and min_pos == snps[i][1]:
                info = snps[i]
                if len(info) <= 2: continue
                #parse out haplotype config info
                ht = map(int, info[9].split('/'))
                #get the configuration
                phased_alleles = [alleles[i][ht[0]], alleles[i][ht[1]]]
                alleles[i] = phased_alleles
        
        #take note that for this position we need to get allele counts in plasma samaples
        loci[min_pos] = alleles
        sp.add_pos(min_pos, data)
        print >>out_pos_file, min_pos, ": M:", alleles[M], " P:", alleles[P], " F:", alleles[F]
            
        #read input: next SNP
        for i in [M, P, F]:
            if min_chr >= snps[i][0] and min_pos >= snps[i][1]:
                snps[i] = in_files[i].readline().split('\t')

        #END WHILE
    
    #set up datastructures for counting allele support in diffrenct SAM files
    posInfo = [dict() for i in ALL]
    for R in [PLASMA, MR, PR]:
        posInfo[R] = copy.deepcopy(data)
        
    #fetch the reads in plasma SAM file and get counts for the positions originally specified in 'data'
    for R in [PLASMA, MR, PR]:
        while True:
            line = in_files[R].readline()
            if not line: break
            if len(line) > 0 and line[0] == '@': continue
            sp.pile_up(sp.mapping_parser(line), posInfo[R])    
    
    #print info / compute stats for each SNP position
    for pos in sorted(data.keys()):
        alleles = loci[pos]
        #if alleles[M][0] != alleles[M][1]:
        #output M, P, F alleles at this SNP locus
        for i in [M, P]:
            print >>out_files[i], alleles[i][0], alleles[i][1]
        print >>out_files[F], alleles[F][0], alleles[F][1], 3
        
        #print the plasma allele counts 
        nuc_counts = posInfo[PLASMA][pos]
        tmp = []
        for nuc in 'ACGT': #to make sure they are in the right order
            try:
                tmp.append(str(nuc_counts[nuc]))
            except KeyError:
                tmp.append('0')
        print >>out_files[PLASMA], ' '.join(tmp)
    
    
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    
    main()

