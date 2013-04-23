#!/usr/bin/python2

import argparse
import subprocess
import random

def main():
    parser = argparse.ArgumentParser(description='Prepare SNP data for FCNV. Read filenames for M, P, F, and plasma .vcf files.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to .vcf files with M, P, F SNPs and plasma mappings')
    args = parser.parse_args()
    
    if len(args.filenames) != 4: die("No enough arguments: missing file names!")
    
    #treat these as CONSTANTS!
    M = 0; P = 1; F = 2; PLASMA = 3;
    ALL = [M, P, F, PLASMA]
    
    #list of input files
    in_files = [open(args.filenames[i], "r" ) for i in ALL]
    
    #list of output files
    out_files = [None for i in ALL]
    out_files[M] = open("M_alleles.txt", "w")
    out_files[P] = open("P_alleles.txt", "w")
    out_files[F] = open("F_alleles.txt", "w")
    out_files[PLASMA] = open("plasma_samples.txt", "w")
    
    #read SNPs from M, P, F vcf files
    snps = [[] for i in ALL]
    for i in [M, P, F]:
        snps[i] = in_files[i].readline().split("\t")
    
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
        
        #for the SNP position, get its alignment info in plasma reads
        while True:
            snps[PLASMA] = in_files[PLASMA].readline().split("\t")
            try:
                pos = int(snps[PLASMA][1])
            except: 
                continue
            if min_pos == pos and min_chr == snps[PLASMA][0]:
                #we found the right position in plasma vcf
                break
            elif min_chr == snps[PLASMA][0]:
                #skip positions we know we don't care about
                for i in xrange(pos - min_pos - 2): plasma_in.readline()
        
        #get alleles
        alleles = [['.', '.'] for x in ALL]
        alleles[PLASMA] = [snps[PLASMA][3], snps[PLASMA][4]]
        for i in [M, P, F]:
            #if there is a SNP in the data at this position, use it
            if min_chr == snps[i][0] and min_pos == snps[i][1]:
                alleles[i] = [snps[i][3], snps[i][4]]

            
        #if there is no info on some allele, impute the reference allele
        ref_allele = alleles[PLASMA][0]
        for i in ALL:
            for a in [0, 1]:
                if alleles[i][a] == '.': alleles[i][a] = ref_allele
        
        #get allele counts in reads
        allele_counts = [[-1, -1] for i in ALL]
        for i in ALL:
            info = snps[i]
            if len(info) <= 2: continue
            #parse out number of supporting reads for reference allele fwd and bck strand, alternate allele fwd and bck strand
            DP = dict( [x.split('=') for x in info[7].split(';')] )['DP4'].split(',')
            DP = map(int, DP)
            allele_counts[i][0] = sum(DP[0:2]) #num of supporting reads for reference allele
            allele_counts[i][1] = sum(DP[2:4]) #num of supporting reads for alternate allele
        
        #if reference allele has no support => the site is homozygous alternative
        for i in ALL:
            if allele_counts[i][0] == 0:
                alleles[i][0] = alleles[i][1]
        
        #ignore homozygous alternative sites that are the same in both M, P, and F:
        ref_support = 0
        for i in [M, P, F]:
            if min_chr == snps[i][0] and min_pos == snps[i][1]:
                ref_support += allele_counts[i][0]
        
        if ref_support != 0:
            #print min_chr, min_pos, "M:", alleles[M], allele_counts[M], "P:", alleles[P], allele_counts[P], \
            # "plasma:", alleles[PLASMA], allele_counts[PLASMA]
            for i in [M, P]:
                print >>out_files[i], alleles[i][0], alleles[i][1]
            print >>out_files[F], alleles[F][0], alleles[F][1], 3
            
            nuc_counts = dict([['A', 0], ['C', 0], ['G', 0], ['T', 0]])
            nuc_counts[alleles[PLASMA][0]] = allele_counts[PLASMA][0]
            nuc_counts[alleles[PLASMA][1]] = allele_counts[PLASMA][1]
            tmp = []
            for nuc in ['A', 'C', 'G', 'T']: #to make sure they are in the right order
                tmp.append(str(nuc_counts[nuc]))
            print >>out_files[PLASMA], " ".join(tmp)
        
        
        #read input: next SNP
        for i in [M, P, F]:
            if min_chr >= snps[i][0] and min_pos >= snps[i][1]:
                snps[i] = in_files[i].readline().split("\t")
    
    
if __name__ == "__main__":
    #import doctest
    #doctest.testmod()
    
    main()

