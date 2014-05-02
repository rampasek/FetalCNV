#!/usr/bin/pypy

import math
import argparse
import sys        

def MS(alleles, j, N, n_plus):
    MSP_j = 0
    MSG_j = 0

    p_j = 0    
    for i in range(N):
        p_j += alleles[i][j]
        
    p_j = float(p_j)/n_plus
    
    for i in range(N):
        n_i = sum(alleles[i])
        p_ij = float(alleles[i][j])/n_i
        MSP_j += n_i*(p_ij - p_j) ** 2
        MSG_j += n_i*p_ij*(1 - p_ij)

    return MSP_j/float(N-1), MSG_j/float(n_plus - N)
    
def rho_mom(alleles):
    
    k = len(alleles[0])

    rho_num = 0
    rho_denom = 0
    N = len(alleles)
    
    n_plus = 0    
    for i in range(N):
        n_plus += sum(alleles[i])    
    
    sum_n_sqr = sum(sum(allele) ** 2 for allele in alleles)
    n_c = (n_plus - sum_n_sqr/float(n_plus)) / float(N-1)
    #print alleles
    
    for j in range(k):
        MSP_j, MSG_j = MS(alleles, j, N, n_plus)
        rho_num += MSP_j - MSG_j
        rho_denom += MSP_j + (n_c - 1) * MSG_j
        #print "\tj {2} MSP {0} MSG {1}".format(MSP_j, MSG_j, j)
    
    #print "rhonum {0} rhodenom {1}".format(rho_num, rho_denom)
    #print k, N
    try:
        rho = rho_num/rho_denom
    except:
        rho = 0

    return rho

def main():
    #parse command line arguments
    parser = argparse.ArgumentParser(description='Compute empirical variance of binomial parameter for allele distribution.')
    parser.add_argument('filenames', type=str, nargs='+', help='files containing binomial parameter stats')
    args = parser.parse_args()
    
    print "#in files: ", len(args.filenames)
    #print args.filenames
    
    parameterStats = dict()
    #return 0
    for fname in args.filenames:
        fin = open(fname, 'r')
        total = int(fin.readline())
        for i in range(total):
            mean = float(fin.readline())
            vals = map(lambda x: tuple(map(int, x.split('/'))), fin.readline().split())
            if mean not in parameterStats: parameterStats[mean] = list()
            parameterStats[mean] += vals
        fin.close()
    
#    #debug output
#    for expP in sorted(parameterStats.keys()):
#        print "----------", expP, "----------"
#        for x in parameterStats[expP]:
#            print x
    sum_val = 0
    num_vals = 0
    #for each category (expected P) compute the MoM estimate of rho
    for expP in sorted(parameterStats.keys()):
        if len(parameterStats[expP]) > 1:
            rho = rho_mom(parameterStats[expP])
            #print "RHO >>>>>", rho
            
            if rho != 0:
                sum_val += abs(rho)*len(parameterStats[expP]) 
                num_vals += len(parameterStats[expP]) 
    
    print "Avg rho:", sum_val / num_vals
    
    
    
if __name__ == "__main__":
    #import doctest
    #doctest.testmod()
    
    main()
    



