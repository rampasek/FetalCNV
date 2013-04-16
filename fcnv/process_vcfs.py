#!/usr/bin/python2

import argparse
import subprocess
import random

def main():
    parser = argparse.ArgumentParser(description='Prepare SNP data for FCNV. Read filenames for M, P, F, and plasma .vcf files.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to .vcf files with M, P, F SNPs and plasma mappings')
    args = parser.parse_args()
    
    if len(args.filenames) != 4: die("No enough arguments: missing file names!")

    M_in = open(args.filenames[0], "r" )
    P_in = open(args.filenames[1], "r" )
    F_in = open(args.filenames[2], "r" )
    plasma_in = open(args.filenames[3], "r" )
    
    M_out = open("M_alleles.txt", "w")
    P_out = open("P_alleles.txt", "w")
    F_out = open("F_alleles.txt", "w")
    plasma_out = open("plasma_samples.txt", "w")
    
    snps = [[], [], []]
    snps[0] = M_in.readline().split("\t")
    snps[1] = P_in.readline().split("\t")
    snps[2] = F_in.readline().split("\t")
    while len(snps[0])>2 or len(snps[1])>2: #while there is a SNP positions in M or P
        #get the position of SNP that occure first [str chromosome, int position]
        for i in [0, 1, 2]:
            if snps[i][0] == '': 
                snps[i][0] = 'chrZ'
                snps[i].append(1e15)
            else:
                snps[i][1] = int(snps[i][1])
        min_chr = min(snps[0][0], snps[1][0])
        min_pos = 1e15
        for i in [0,1]:
            if min_chr == snps[i][0] and snps[i][1] < min_pos:
                min_pos = snps[i][1]
        
        #for the SNP position, get its alignment info in plasma reads
        while True:
            plasma_pos_info = plasma_in.readline().split("\t")
            if min_chr == plasma_pos_info[0] and min_pos == int(plasma_pos_info[1]): break
        
        #get alleles
        M_alleles = ['.', '.']
        P_alleles = ['.', '.']
        F_alleles = ['.', '.']
        plasma_alleles = [plasma_pos_info[3], plasma_pos_info[4]]
        if min_chr == snps[0][0] and min_pos == snps[0][1]:
            M_alleles = [snps[0][3], snps[0][4]]
        if min_chr == snps[1][0] and min_pos == snps[1][1]:
            P_alleles = [snps[1][3], snps[1][4]]
        if min_chr == snps[2][0] and min_pos == snps[2][1]:
            F_alleles = [snps[2][3], snps[2][4]]
            
        #if there is no info on some allele, impute the reference allele
        ref_allele = plasma_alleles[0]
        for i in [0,1]:
            if M_alleles[i] == '.': M_alleles[i] = ref_allele
            if P_alleles[i] == '.': P_alleles[i] = ref_allele
            if F_alleles[i] == '.': F_alleles[i] = ref_allele
        
        #get allele counts in reads
        allele_counts = [[-1, -1] for x in range(3)]
        for i, info in enumerate([snps[0], snps[1], plasma_pos_info]):
            if len(info) <= 2: continue
            DP = dict([x.split('=') for x in info[7].split(';')])['DP4'].split(',')
            DP = map(int, DP)
            allele_counts[i][0] = sum(DP[0:2])
            allele_counts[i][1] = sum(DP[2:4])
        
        #ignore homozygous alternative sites that are the same in both M, P, and F:
        ref_support = allele_counts[2][0]
        for i in range(2):
            if min_chr == snps[i][0] and min_pos == snps[i][1]:
                ref_support += allele_counts[i][0]
        
        if ref_support != 0:
            #print min_chr, min_pos, "M:", M_alleles, allele_counts[0], "P:", P_alleles, allele_counts[1], "plasma:", plasma_alleles, allele_counts[2]
            print >>M_out, M_alleles[0], M_alleles[1]
            print >>P_out, P_alleles[0], P_alleles[1]
            print >>F_out, F_alleles[0], F_alleles[1], 3
            nuc_counts = dict([['A', 0], ['C', 0], ['G', 0], ['T', 0]])
            nuc_counts[plasma_alleles[0]] = allele_counts[2][0]
            nuc_counts[plasma_alleles[1]] = allele_counts[2][1]
            tmp = []
            for nuc in ['A', 'C', 'G', 'T']: #to make sure they are in the right order
                tmp.append(str(nuc_counts[nuc]))
            print >>plasma_out, " ".join(tmp)
        
        
        #read input: next SNP
        for i in [0, 1, 2]:
            if min_chr >= snps[i][0] and min_pos >= snps[i][1]:
                if i==0: snps[i] = M_in.readline().split("\t")
                if i==1: snps[i] = P_in.readline().split("\t")
                if i==2: snps[i] = F_in.readline().split("\t")
    
    
if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    
    main()

