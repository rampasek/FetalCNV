#!/usr/bin/pypy

import math
import argparse
import sys        

def readParams(param_file):
    """Read the CRF model parameters from file
    
    arguments:
    param_file -- file handler opened for read
    returns dictionary with parameters    
    """
    crfParams = {}
    if not param_file:
        print "readParams: Error, file handler not valid"
        return crfParams
    
    for line in param_file.readlines():
        line = line.rstrip('\n').split('=')
        if len(line) != 2: continue
        line[0] = line[0].strip().rstrip()
        
        if line[0] == 'unaryWeights':
            unaryWeights = map(float, line[1].split())
            crfParams['unaryWeights'] = unaryWeights
            
        elif line[0] == 'binaryWeights':
            binaryWeights = map(float, line[1].split()) #map(lambda x: math.log(float(x)), line[1].split())
            crfParams['binaryWeights'] = binaryWeights
            
        elif line[0] == 'epsWeight':
            epsWeight = float(line[1])
            crfParams['epsWeight'] = epsWeight
            
        elif line[0] == 'sigmaSqr':
            sigmaSqr = float(line[1])
            crfParams['sigmaSqr'] = sigmaSqr
            
        elif line[0] == 'omega':
            omega = float(line[1])
            crfParams['omega'] = omega
            
        else:
            print "Unknown parameter keyword", line[0], line[1]
    
    if len(crfParams) != 5:
        print "readParams: Not enough parameters!!"
    return crfParams

def main():
    #parse command line arguments
    parser = argparse.ArgumentParser(description='Compute average of the input sets of CRF weight parameters.')
    parser.add_argument('outParamFile', type=str, nargs=1, help='file name for the resulting average CRF parameters')
    parser.add_argument('filenames', type=str, nargs='+', help='files containing CRF weights')
    args = parser.parse_args()
    
    out_file_name = args.outParamFile[0]
    print "#in files: ", len(args.filenames)
    if len(args.filenames) < 1:
        print "ERROR: No input CRF parameter files"
        return -1
    #print args.filenames
    
    allParams = list()
    #return 0
    for fname in args.filenames:
        fin = open(fname, 'r')
        allParams.append(readParams(fin))
    
    N = float(len(allParams))
    params = allParams[0]
    for p in allParams[1:]:
        for i in len(p['unaryWeights']):
            params['unaryWeights'][i] += p['unaryWeights'][i]
        for i in len(p['binaryWeights']):
            params['binaryWeights'][i] += p['binaryWeights'][i]
    
    #average
    for i in len(params['unaryWeights']):
        params['unaryWeights'][i] /= N
    for i in len(params['binaryWeights']):
        params['binaryWeights'][i] /= N
    
    #save the trained parameters to the file
    res_param_file = open(out_file_name, "w")    
    for p in sorted(params):
        if isinstance(params[p], list):
            print >>res_param_file, p, "=", " ".join(map(str, params[p]))
        else:
            print >>res_param_file, p, "=", params[p]
    res_param_file.close()
    
    
if __name__ == "__main__":
    #import doctest
    #doctest.testmod()
    
    main()
    



