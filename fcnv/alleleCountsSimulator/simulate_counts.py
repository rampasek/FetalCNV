#!/usr/bin/python2

import math
import argparse
import numpy as np
import sys
import os
import itertools

allele_file_header = ''

def readInput(in_file_name):
    """Read the input from files.
    Read maternal alleles, paternal alleles, and samples per position.
    
    Arguments:
    file name -- string
    returns lists of parsed data    
    
    """
    global allele_file_header
    in_file = open(in_file_name, 'r')
    positions = []
    samples = []
    M = []; P = [];
    MC = []; PC = [];
    while True:
        line = in_file.readline()
        if not line: break
        if line[0] == '#':
            allele_file_header += line
            continue  #skip comment
        line = line.rstrip('\n').split('\t')
        
        #genomic positions and allele support in plasma samples
        positions.append(int(line[0]))
        samples.append(tuple(map(int, line[1:5])))
        
        #maternal and paternal alleles
        M.append(tuple(line[5:7]))
        MC.append(tuple(map(float, line[7:9])))
        
        P.append(tuple(line[9:11]))
        PC.append(tuple(map(float, line[11:13])))     
    
    in_file.close()
    return positions, samples, M, P, MC, PC

def writeOutput(out_file_name, positions, samples, M, P, MC, PC):
    length_ok = len(positions) == len(samples) == len(M) == len(P) == len(MC) == len(PC)
    if length_ok == False:
        print 'Error: output list are not of the same length'
        return
    
    global allele_file_header
    out_file = open(out_file_name, 'w')
    print >>out_file, allele_file_header.rstrip()
    for i in range(len(positions)):
        line = ''
        #genomic positions and allele support in plasma samples
        line += str(positions[i])
        line += '\t{0}\t{1}\t{2}\t{3}'.format(samples[i][0], samples[i][1], samples[i][2], samples[i][3])
        
        #maternal and paternal alleles
        line += '\t{0}\t{1}'.format(M[i][0], M[i][1])
        line += '\t{0}\t{1}'.format(MC[i][0], MC[i][1])
        
        line += '\t{0}\t{1}'.format(P[i][0], P[i][1])
        line += '\t{0}\t{1}'.format(PC[i][0], PC[i][1])  
        
        print >>out_file, line
    out_file.close()

def generatePPstates():
    inheritance_patterns = []
    for mats in range(3):
        for pats in range(3):
            if mats+pats in [0,4]: continue
            inheritance_patterns.append((mats, pats))
            
    #generate phased variants of the inheritance patterns and list of HMM States
    states = []
    phased_patterns = []
    for ip in inheritance_patterns:
        for mPhased in itertools.combinations_with_replacement([0,1], ip[0]):
            for pPhased in itertools.combinations_with_replacement([0,1], ip[1]):
                phased_pattern = (tuple(mPhased), tuple(pPhased))
                phased_patterns.append(phased_pattern)
                #new_state = HMMState(ip, phased_pattern)
                #states.append(new_state)
    #for x in phased_patterns:
    #    print x
    return phased_patterns


def allelesDistributionExtended(Malleles, Mcounts, Palleles, Pcounts, pattern, mix):
    """Compute nucleotides distribution in maternal plasma for a position with
    given maternal and fetal alleles, assuming the given fetal admixture.
    
    Arguments:
    maternal -- maternal alleles
    paternal -- paternal alleles
    mix -- fetal admixture ratio
    returns list(floats) -- list of nucleotides probabilities
    """
    
    nucleotides = ['A', 'C', 'G', 'T']
    
    numFalleles = float(len(pattern[0]) + len(pattern[1]))
    numMinherited = len(pattern[0])
    numPinherited = len(pattern[1])
    numMalleles = float(len(Malleles))
    sumMcount = float(sum(Mcounts))
    sumPcount = float(sum(Pcounts))
    
    #adjust mixture to currently considered event (deletion, duplication, normal)
    adjusted_fetal_admix = mix/2. * numFalleles
    adjusted_maternal_admix = (1.-mix)/2. * numMalleles
    cmix = adjusted_fetal_admix / (adjusted_maternal_admix + adjusted_fetal_admix)
    
    dist = {}
    for nuc in nucleotides: dist[nuc] = 0.
    alphaM = 32/2. + sum(Mcounts)/2.
    alphaP = 39/2. + sum(Pcounts)/2.
    
    #fraction that is trully from mother DNA
    for i, nuc in enumerate(Malleles):
        dist[nuc] += (alphaM + Mcounts[i]) / float(2*alphaM + sumMcount) * (1.-cmix)
        
    #fraction that is fetal but maternaly inherited
    for mHpatt in pattern[0]:
        nuc = Malleles[mHpatt]
        if numMinherited == 2 and pattern[0][0] != pattern[0][1]:
            #if both haplotypes are inherited, use the maternal seq. ratio
            dist[nuc] += (alphaM + Mcounts[mHpatt]) / float(2*alphaM + sumMcount) * (cmix*2./numFalleles)
        else:
            dist[nuc] += cmix / numFalleles
        
    #fraction that is fetal but paternaly inherited
    for pHpatt in pattern[1]:
        nuc = Palleles[pHpatt]
        dist[nuc] += cmix / numFalleles
    
    dist_list = [dist[nuc] for nuc in nucleotides]
    return dist_list

def simulateCounts(samples, M, P, gt, ff):
    newS = []
    phased_patterns = generatePPstates()
    
    for i in range(len(samples)):
        ps = allelesDistributionExtended(M[i], (30,30), P[i], (30,30), phased_patterns[gt[i]], ff)
        sim = np.random.multinomial(sum(samples[i]), ps)
        newS.append(sim)
        
        #print samples[i], M[i], P[i], ps
        #print tuple(sim), '\n'
        #if i == 5: break
    
    return newS
        
def main():
    #parse command line arguments
    parser = argparse.ArgumentParser(description='Performs fetal CNV analysis from maternal plasma and phased parental data.')
    parser.add_argument('input', type=str, nargs=1, help='path to input file with allele counts in plasma and parental haplotypes')
    parser.add_argument('target', type=str, nargs=1, help='path to file with background truth - "target"')
    parser.add_argument('outDir', type=str, nargs=1, help='path to output directory')
    parser.add_argument('--ff', type=float, help='fetal mixture ratio')
    args = parser.parse_args()
    
    in_file_name = args.input[0]
    target_file_name = args.target[0]
    out_dir = args.outDir[0]
    if args.ff == None:
        print "ERROR: no fetal fraction specified!!!"
        return -1
    ff = args.ff
    
    #print input info
    print "------------------------------------------"
    print "Running simulator, input parameters:"
    print "input:", in_file_name
    print "target:", target_file_name
    print "outDir:", out_dir
    print "--ff:", ff
    print "------------------------------------------"
    os.system("hostname")
    
    #read the pre-processed input
    snp_positions, samples, M, P, MSC, PSC = readInput(in_file_name)
    
    #snp_positions = []
    ground_truth = []
    target_file = open(target_file_name, 'r')
    for line in target_file.readlines():
        line = line.rstrip("\n").split("\t")
        #snp_positions.append(int(line[0]))
        ground_truth.append(int(line[-1]))
    target_file.close()
    
    
    file_name_prefix = target_file_name.split('/')[-1].split('.')[0].replace(':', '-')
    out_file_name = out_dir + '/' + file_name_prefix + '.sim_allele_counts.txt'
    #print out_file_name
    #out_file = open(out_file_name, 'w')
    
    #print snp_positions[0], samples[0], M[0], P[0], MSC[0], PSC[0], ground_truth[0]
    
    MSC = [(30,30)] * len(MSC)
    PSC = [(30,30)] * len(MSC)
    samplesSim = simulateCounts(samples, M, P, ground_truth, ff)
    
    writeOutput(out_file_name, snp_positions, samplesSim, M, P, MSC, PSC)
    
if __name__ == "__main__":
    #import doctest
    #doctest.testmod()
    
    main()
    