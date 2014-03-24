#!/usr/bin/pypy

import argparse
import sys
import random
import copy
from datetime import datetime
import pplabelsHMM

def parseBeagleVCF(infile, ch):
    loci = dict()
    processed_chr = ''
    for line in infile.readlines():
        line = line.rstrip('\n')
        if line == '': continue
        #skip the header
        if len(line) > 0 and line[0] == '#': continue
        #split
        line = line.split('\t')
        
        snp = [0] * 5
        #get chromosome
        if processed_chr == '': processed_chr = line[0]
        if processed_chr != line[0]: print "WARNING: multiple chromosomes in the VCF", processed_chr, "|", line[0]
        #position
        snp[0] = int(line[1])
        #dbSNP rsid
        #snp[x] = line[2]
        #ref and alt alleles
        snp[1] = (line[3], line[4])

        if ch == True:
            #parse out CHILD haplotype config info
            ht = map(int, line[11].split(':')[0].split('|'))
            #hapA and hapB
            snp[2] = ht #(snp[2][ht[0]], snp[2][ht[1]])
        else:
            #parse out MATERNAL haplotype config info
            ht = map(int, line[9].split(':')[0].split('|'))
            #hapA and hapB
            snp[2] = ht #(snp[2][ht[0]], snp[2][ht[1]])
            
            #parse out PATERNAL haplotype config info
            ht = map(int, line[10].split(':')[0].split('|'))
            #hapA and hapB
            snp[3] = ht #(snp[2][ht[0]], snp[2][ht[1]])
        
        loci[snp[0]] = tuple(snp)
    infile.close()
    return loci


def main():
    #parse command line arguments
    parser = argparse.ArgumentParser(description='Annotate the children haplotypes by inheritance origin (M_A, M_B, P_A, P_B) -- 8 possible states')
    parser.add_argument('trioPhase', type=str, nargs=1, help='path to Beagle VCF file with pedigree phasing')
    parser.add_argument('mpPhase', type=str, nargs=1, help='path to Beagle VCF file with maternal and paternal haplotypes')
    parser.add_argument('outFile', type=str, nargs=1, help='name of output file for haplotype annotations')
    args = parser.parse_args()
    
    trio_phase_file = open(args.trioPhase[0], "r")
    mp_phase_file = open(args.mpPhase[0], "r")
    
    ### alleles and haplotypes for particular SNP positions
    #parse Beagle VCF file from TRIO phasing
    lociCH = parseBeagleVCF(trio_phase_file, True)
    #parse Beagle VCF file from maternal and paternal phasing
    
    #mp_phase_file = open(args.trioPhase[0], "r") #to test, use the pedigree M, P phasing
    lociMP = parseBeagleVCF(mp_phase_file, False)

    #for i in range(25):
    #    print lociCH.items()[i], '\n', lociMP.items()[i], '\n'
    
    #delete SNPs from child's haplotypes that are not in MP phasing
    ch_pos2del = []
    for pos in lociCH.keys():
        if pos not in lociMP:
            ch_pos2del.append(pos)
    for pos in ch_pos2del:
        del lociCH[pos]
    
    #insert dummy SNPs to child's haplotypes at positions that are MP phasing but not in CH phasing
    for pos in lociMP.keys():
        if pos not in lociCH:
            lociCH[pos] = tuple([-1, tuple(['x','x']), tuple([-1,-1])])
            
    if len(lociCH.keys()) != len(lociMP.keys()):
        print "ERROR!!!!!!!!!! list length do not match"
        return 0
        
    keysCH = list(sorted(lociCH.keys()))
    #keysMP = list(sorted(lociMP.keys()))
    
    #print len(keysCH), len(keysMP)
    #for i in range(len(keysCH)):
    #    if keysCH[i] != keysMP[i]: print "pos", i
    #    if lociCH[keysCH[i]][1] != lociMP[keysMP[i]][1]: print "loci err:", lociCH[keysCH[i]], lociMP[keysMP[i]]
        
    ch_haps = [None] * len(keysCH)
    mp_haps = [None] * len(keysCH)
    for i in range(len(keysCH)):
        ch_haps[i] = lociCH[keysCH[i]][2]
        mp_haps[i] = tuple(lociMP[keysCH[i]][2:4])
    
    #for i in range(200):
    #    print ch_haps[i], mp_haps[i]
    
    
    hmm = pplabelsHMM.FCNV(keysCH)
    vp, vp_states = hmm.viterbiPath(ch_haps, mp_haps)
    
    #snp positions and their haloptype annotation
    out_file = open(args.outFile[0], "w")
    for i in range(len(keysCH)):
        print >>out_file, keysCH[i], vp[i]
    out_file.close()
    
    ## DEBUG output and stats
#    print len(keysCH), len(vp)
#    
#    switches = 0
#    leng = 0
#    prev = -1
#    for x in vp:
#        leng += 1
#        if x!= prev:
#            switches += 1
#            print prev, '(', leng, ')', ' ->', x
#            leng = 0
#        prev = x
#    print prev, '(', leng, ')'
#    print "switches:", switches
#    
#    #print vp_states
#    #print '\n'.join(map(str, vp_states))
    
    
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    
    main()

