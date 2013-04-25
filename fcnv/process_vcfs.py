#!/usr/bin/python2

import argparse
import subprocess
import random

#SAM bitewise flags
ALL_PROPERLY_ALIGNED = 0x2
IS_UNMAPPED = 0x4
MP_UNMAPPED = 0x8
#CIGAR operators
cigar_qr_op = 'HSIM=X'
cigar_db_op = 'NDM=X'
cigar_clip_op = 'HS'

def mapping_parser(m):
    '''
    Parse a read in SAM format, return a dictionary with filled in fields of interest.
    '''
    if isinstance(m, str):
        m = m.strip().split('\t')
    d = {}
    d['flag'] = int(m[1])   # flags
    d['chr'] = m[2]         # chr
    d['pos'] = int(m[3])    # pos
    d['cigar'] = m[5]       # cigar string
    d['seq'] = m[9]         # sequence
    
    return d

def parse_cigar_string(s):
    '''
    Parse given CIGAR string to a list of operators.
    '''
    res = []
    crt_len = ''
    i = 0
    while i < len(s):
        if str.isdigit(s[i]):
            crt_len += s[i]
        else:
            res.append([s[i], int(crt_len)])
            crt_len = ''
        i += 1
    return res
    
def pile_up(read, data):
    '''
    Store allele evidence given by @read.
    '''
    #take only reads mapped in proper pair
    if read['flag'] & ALL_PROPERLY_ALIGNED == 0: return
    
    pos_qr = 0
    pos_db = read['pos']
    op = parse_cigar_string(read['cigar'])
    for o in op:
        if o[0] == 'H': continue
        elif o[0] in 'SI': pos_qr += o[1]
        elif o[0] in 'ND': pos_db += o[1]
        elif o[0] in 'M=X':
            for i in range(o[1]):
                add_support(pos_db, read['seq'][pos_qr], data)
                pos_db += 1
                pos_qr += 1

def add_pos(pos, data):
    '''
    Add key @pos to @data.
    '''
    pos = str(pos)
    if not pos in data: data[pos] = dict()

def add_support(pos, nuc, data):
    '''
    Increase counter for nucleotide @nuc at position @pos, motify datastructure @data.
    Ignore positions @pos that are not in @data.
    '''
    pos = str(pos)
    if not pos in data: return #if we don't need info for this position
    if not nuc in 'ACGT': return
    try:
        data[pos][nuc] += 1
    except:
        data[pos].update({nuc:1})

def main():
    #parse ARGs
    parser = argparse.ArgumentParser(description='Prepare SNP data for FCNV. Read filenames for M, P, and F .vcf files and plasma .sam file.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to .vcf files with M, P, F SNPs and *sorted* plasma reads in SAM format')
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
    
    #allele counts in plasma samples for particular positions
    data = dict()
    
    #read SNPs from M, P, F vcf files
    snps = [[] for i in [M, P, F]]
    for i in [M, P, F]:
        snps[i] = in_files[i].readline().split('\t')
    
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
        
        #get allele counts in reads
        allele_counts = [[-1, -1] for i in [M, P, F]]
        for i in [M, P, F]:
            info = snps[i]
            if len(info) <= 2: continue
            #parse out number of supporting reads for reference allele fwd and bck strand, alternate allele fwd and bck strand
            DP = dict( [x.split('=') for x in info[7].split(';')] )['DP4'].split(',')
            DP = map(int, DP)
            allele_counts[i][0] = sum(DP[0:2]) #num of supporting reads for reference allele
            allele_counts[i][1] = sum(DP[2:4]) #num of supporting reads for alternate allele
        
        #if reference allele has no support => the site is homozygous alternative
        for i in [M, P, F]:
            if allele_counts[i][0] == 0:
                alleles[i][0] = alleles[i][1]
        
        #ignore homozygous alternative sites that are the same in both M, P, and F:
        #ref_support = 0
        #for i in [M, P, F]:
        #    if min_chr == snps[i][0] and min_pos == snps[i][1]:
        #        ref_support += allele_counts[i][0] 
        #if ref_support != 0:
            #print min_chr, min_pos, "M:", alleles[M], allele_counts[M], "P:", alleles[P], allele_counts[P], \
            # "plasma:", alleles[PLASMA], allele_counts[PLASMA]
            
        #output M, P, F alleles at this SNP locus
        for i in [M, P]:
            print >>out_files[i], alleles[i][0], alleles[i][1]
        print >>out_files[F], alleles[F][0], alleles[F][1], 3
        
        #take note that for this position we need to get allele counts in plasma samaples
        add_pos(min_pos, data)
        
        #read input: next SNP
        for i in [M, P, F]:
            if min_chr >= snps[i][0] and min_pos >= snps[i][1]:
                snps[i] = in_files[i].readline().split('\t')
    
    #read the reads in plasma SAM file and get counts for the positions specified in 'data'
    while True:
        line = in_files[PLASMA].readline()
        if not line: break
        pile_up(mapping_parser(line), data)
    
    #print the plasma allele counts    
    for pos in sorted(data.keys()):
        nuc_counts = data[pos]
        tmp = []
        for nuc in ['A', 'C', 'G', 'T']: #to make sure they are in the right order
            try:
                tmp.append(str(nuc_counts[nuc]))
            except KeyError:
                tmp.append('0')
        print >>out_files[PLASMA], ' '.join(tmp)
    
    
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    
    main()

