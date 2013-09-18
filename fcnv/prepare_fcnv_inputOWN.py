#!/usr/bin/python2
#example usage: time prepare_fcnv_input.py mp.phase.vcf __plasma.dup100kb.sam __M.part.sam __P.part.sam /dupa-filer/laci/centromeres >log_prepare_fcnv_input.txt 2>>log_prepare_fcnv_input.txt &
import argparse
from sys import exit
import random
import copy
import samParse as sp
from datetime import datetime

def is_within_intervals(num, intervals):
    for interval in intervals:
        if interval[0] <= num <= interval[1]:
            return True
    return False

def main():
    #parse ARGs
    parser = argparse.ArgumentParser(description='Prepare SNP support data for FCNV. Read filenames: for joined M&P phased .vcf file; plasma, M, and P .sam files; and for centromeres list.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to 1) .vcf file with phased M & P SNPs; 2) reads in SAM format for plasma, M, and P samples; 3) centromeres list file.')
    args = parser.parse_args()
    
    if len(args.filenames) != 5: exit("Unexpected number of arguments passed! Expecting 5 filenames.")
    
    #treat these as CONSTANTS!
    MP = 0; PLR = 1; MR = 2; PR = 3; CT = 4; #in_files
    ALL = [MP, PLR, MR, PR, CT]
    M = 0; P = 1; #maternal, paternal
    ALDOC = 0; GT = 1; #out_files: allele DOC and ground truth
    
    #list of input files
    in_files = [open(args.filenames[i], "r" ) for i in ALL]
    
    #read centromeres positions
    centromeres = dict()
    for line in in_files[CT].readlines():
        line = line.rstrip('\n').split('\t')
        if line[0] not in centromeres.keys(): centromeres[line[0]] = []
        centromeres[line[0]] += [(int(line[1]), int(line[2]))]
        
    
    #allele counts in plasma samples for particular positions
    pos_data = dict()
    loci = dict()
    processed_chr = ''
    skipped_in_centromere = 0
    
    print "  Getting union of SNP positions and corresponding list of alleles " + datetime.now().strftime('%m-%d-%H-%M')
    #read SNPs from M, P, F vcf files
    snps = [[] for i in [M, P]]
    
    #get genotypes for all positions in UNION of M and P SNP positions
    while True:
        line = in_files[MP].readline()
        #skip if part of the header
        while len(line) > 0 and line[0] == '#': line = in_files[MP].readline()
        if not line: break
        
        fields = line.rstrip('\n').split('\t')
        if processed_chr == '': processed_chr = 'chr'+fields[0]
        if processed_chr != 'chr'+fields[0]:
            print "WARNING: multiple chromosomes in the input", processed_chr, "|", 'chr'+fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]
        
        #get M and P haplotypes
        snps[M] = fields[9].split(':')[0].split('|')
        for x in [0, 1]:
            if snps[M][x] == '0': snps[M][x] = ref
            else: snps[M][x] = alt
            
        snps[P] = fields[10].split(':')[0].split('|')
        for x in [0, 1]:
            if snps[P][x] == '0': snps[P][x] = ref
            else: snps[P][x] = alt
        
        #if in centromere region, skip
        centromere_regions = centromeres[processed_chr]
        if is_within_intervals(pos, centromere_regions):
            skipped_in_centromere += 1
            continue
        
        #take note that for this position we need to get allele counts in plasma samaples
        alleles = (snps[M], snps[P])
        loci[pos] = alleles
        sp.add_pos(pos, pos_data)
        
        #END WHILE
    
    print "  Piling up the reads " + datetime.now().strftime('%m-%d-%H-%M')
    #set up datastructures for counting allele support in diffrenct SAM files
    posInfo = [dict() for i in ALL]
    for R in [PLR, MR, PR]:
        posInfo[R] = copy.deepcopy(pos_data)
        
    #fetch the reads in plasma SAM file and get counts for the positions originally specified in 'pos_data'
    for R in [PLR, MR, PR]:
        while True:
            line = in_files[R].readline()
            if not line: break
            if len(line) > 0 and line[0] == '@': continue
            sp.pile_up(sp.mapping_parser(line), posInfo[R])    
    
    
    print "  Writing output " + datetime.now().strftime('%m-%d-%H-%M')
    #list of output files
    out_files = [None for i in [ALDOC, GT]]
    out_files[ALDOC] = open(processed_chr + "_alleles_docOWN.txt", "w")
    out_files[GT] = open(processed_chr + "_targetOWN.txt", "w")
    print >>out_files[ALDOC], '#POS\tA\tC\tG\tT\tM_hapA\tM_hapB\tDP_hapA\tDP_hapB\tP_hapA\tP_hapB\tDP_hapA\tDP_hapB'
    
    skipped_low_doc = 0
    #print info / compute stats for each SNP position
    for pos in sorted(pos_data.keys()):
        alleles = loci[pos]
        
        #print the plasma allele counts 
        nuc_counts = posInfo[PLR][pos]
        tmp = []
        for nuc in 'ACGT': #to make sure they are in the right order
            try:
                tmp.append(str(nuc_counts[nuc]))
            except KeyError:
                tmp.append('0')
                
        #if the plasma coverage is too low, skip this position        
        if sum(map(int, tmp)) < 50: 
            print pos, "- low overall coverage", sum(map(int, tmp))
            skipped_low_doc += 1
            continue
        
        print >>out_files[ALDOC], str(pos) + '\t' + '\t'.join(tmp),
        
        #output M, P alleles at this SNP locus
        for i, r in [(M, MR), (P, PR)]:
            a1 = alleles[i][0]
            a2 = alleles[i][1]
            count_a1 = 0
            count_a2 = 0
            try: count_a1 = posInfo[r][pos][a1]
            except: print i, pos, a1, posInfo[r][pos], alleles[i]
            try: count_a2 = posInfo[r][pos][a2]
            except: print i, pos, a2, posInfo[r][pos], alleles[i]
            
            if a1 == a2:
                count_a1 /= 2.
                count_a2 /= 2.
            
            print >>out_files[ALDOC], '\t{0}\t{1}\t{2}\t{3}'.format(a1, a2, count_a1, count_a2),
        print >>out_files[ALDOC], '\n',
        
        if 10181440 <= pos and pos <= 10281440:
            print >>out_files[GT], '{0}\t{1}\t{2}\t{3}'.format(pos, 'N', 'N', 6)
        else:
            print >>out_files[GT], '{0}\t{1}\t{2}\t{3}'.format(pos, 'N', 'N', 3)
     
    print "Low overall coverage positions ignored:", skipped_low_doc
    print "Ignored positions in centromere regions:", skipped_in_centromere   
    print "DONE " + datetime.now().strftime('%m-%d-%H-%M')
    
if __name__ == '__main__':
    main()

