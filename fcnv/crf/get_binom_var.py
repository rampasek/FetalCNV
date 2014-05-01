#!/usr/bin/pypy

import math
import argparse
import sys        

def MSP(j):
    pass
    
def MSG(j):
    pass

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
    
    #debug output
    for expP in sorted(parameterStats.keys()):
        print "----------", expP, "----------"
        for x in parameterStats[expP]:
            print x
    
    #for each category (expected P) compute the MoM estimate of rho
    for expP in sorted(parameterStats.keys()):
        pass
    
    
#    parameterStats = dict()
#    #return 0
#    for fname in args.filenames:
#        fin = open(fname, 'r')
#        total = int(fin.readline())
#        for i in range(total):
#            mean = float(fin.readline())
#            vals = map(float, fin.readline().split())
#            if mean not in parameterStats: parameterStats[mean] = list()
#            parameterStats[mean] += vals
#        fin.close()
#    
#    summ = 0.
#    total = 0.
#    catBegin = -1
#    
##    for expP in sorted(parameterStats.keys()):
##        print "----------", expP, "----------"
##        
##        for x in parameterStats[expP]:
##            v = (x-expP)**2v
##            N = 78
##            rho = (N*v)/(1. + (N-1)) * 1./(expP*(1. - expP))
##            print rho,
##        print "\n================="
#    
#    for expP in sorted(parameterStats.keys()):
#        #empAvg = sum(parameterStats[expP]) / float(len(parameterStats[expP]))
#        #empVar = sum([ (x-empAvg)**2 for x in parameterStats[expP]]) / float(len(parameterStats[expP]))
#        #expVar = sum([ (x-expP)**2 for x in parameterStats[expP]]) / float(len(parameterStats[expP]))
#        #print expP, ':  ', empAvg, empVar, math.sqrt(empVar), " against exp.:", expVar, math.sqrt(expVar)
#        
#        if catBegin == -1: catBegin = expP
#        if expP - catBegin > 0.03:
#            #print '[', catBegin, expP,') :', summ/total, '      total:', total
#            if total != 0:
#                print 'if {0} <= x < {1}: var = {2}'.format(catBegin, expP, summ/total)
#                catBegin = expP
#                summ = 0.
#                total = 0.
#        
#        #compute the variance 
#        leng = len(parameterStats[expP])
#        vals = sorted(parameterStats[expP])[leng/10:9*leng/10]
#        
#        summ += sum([ (x-expP)**2 for x in vals])
#        total += len(vals)
#    
#    #print summ/total, '      total:', total
#    print 'if {0} <= x < {1}: var = {2}'.format(catBegin, expP, summ/total)
    
if __name__ == "__main__":
    #import doctest
    #doctest.testmod()
    
    main()
    



