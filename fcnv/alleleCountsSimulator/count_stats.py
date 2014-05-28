#!/usr/bin/env python2.7

import math
import argparse
import numpy as np
import sys
import os
import itertools
from scipy.interpolate import UnivariateSpline
from matplotlib import pyplot as plt

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

def computeCountsStats(samples, M, P, MSC, PSC, gt, ff):
    stats = [ dict() for i in range(30) ] 
    nucleotides = ['A', 'C', 'G', 'T']
    phased_patterns = generatePPstates()
    
    for sid in range(len(samples)):
        ps = allelesDistributionExtended(M[sid], (30, 30), P[sid], (30, 30), phased_patterns[gt[sid]], ff)
        nuc_counts = samples[sid]
        
        expinherited = [0] * 4
        for x in M[sid] + P[sid]:
            expinherited[nucleotides.index(x)] = 1
        
        seen = [0] * 4
        for x in range(4):
            if nuc_counts[x] != 0: seen[x] = 1
            
        for x in range(4):
            seen[x] = int(seen[x] or expinherited[x])
        
        als = [i for i in range(4) if seen[i] == 1]
        
        #print nuc_counts, '|', als, '|', seen
        if len(als) > 2: continue
        if len(als) == 1:
            for x in range(4):
                if x != als[0]: 
                    als.append(x)
                    break
        
        m = len(als)
        eps = 0.005 #.5/sum(nuc_counts)
        dmps = [0.] * m
        addEps = False
        for i in range(m):
            if ps[als[i]] < 0.0001: addEps = True
        for i in range(m):
            dmps[i] = ps[als[i]] + eps*int(addEps)
        dmps = [ x/sum(dmps) for x in dmps]
        
        #make the one with higher p to be the first one
        if dmps[1] > dmps[0]:
            dmps[0], dmps[1] = dmps[1], dmps[0]
            als[0], als[1] = als[1], als[0]
        
        nc = [nuc_counts[als[i]] for i in range(m)]
        #print nc, dmps, " | ", nuc_counts, ps, list(np.array(nuc_counts)/float(sum(nuc_counts)))
        
        if tuple(dmps) not in stats[gt[sid]]: stats[gt[sid]][tuple(dmps)] = list()
        stats[gt[sid]][tuple(dmps)].append(tuple(nc))
    
    return stats
        
def main():
    #parse command line arguments
    parser = argparse.ArgumentParser(description='')
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

    stats = computeCountsStats(samples, M, P, MSC, PSC, ground_truth, ff)
    #print stats
    
    num = 0
    phased_patterns = generatePPstates()
    for state in range(30):
        for mu, emp in stats[state].iteritems():
            print state, phased_patterns[state], ":", mu, len(emp)
            #continue
            plt.figure(figsize=(16, 12))
            n = 150
            bins = np.linspace(0., 1., n)
            
            s = [float(e[0])/sum(e) for e in emp]
            plt.hist(s, bins, color='blue', alpha=0.5, label='empirical')
            
            #p, x = np.histogram(s, bins=n) # bin it into n bins
            #x = x[:-1] + (x[1] - x[0])/2   # convert bin edges to centers            
            #f = UnivariateSpline(x, p, s=1)
            #plt.plot(x, f(x))
            
            s = [np.random.multinomial(sum(e), mu) for e in emp]
            s = [float(e[0])/sum(e) for e in s]
            plt.hist(s, bins, color='green', alpha=0.5, label='sinthetic')
            
            #p, x = np.histogram(s, bins=n) # bin it into n bins
            #x = x[:-1] + (x[1] - x[0])/2   # convert bin edges to centers
            #f = UnivariateSpline(x, p, s=1)
            #plt.plot(x, f(x))
            
            plt.legend(loc='upper left')
            plt.axvline(mu[0], color='r', linestyle='dashed', linewidth=2)
            
            #plt.show()
            num += 1
            plt.savefig('hist/'+str(num)+'hist'+str(state)+str(phased_patterns[state])+'.png', dpi = 100)
            plt.close('all')
    
    #writeOutput(out_file_name, snp_positions, samplesSim, M, P, MSC, PSC)
    
if __name__ == "__main__":
    #import doctest
    #doctest.testmod()
    
    main()
    
