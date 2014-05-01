#!/usr/bin/pypy

#time pypy fcnv2.py IM1-A-chr1-29032136-30032136-duplicate.alleles_doc.txt IM1-A-chr1-29032136-30032136-duplicate.pplabels.txt /dev/null /dev/null /dev/null crfParams.in --trainGrad crfParams.outMM  > logfcnv2.log 2>&1 &

import math
import argparse
#import numpypy as np
import sys
import fcnvCRF
from datetime import datetime
import os


def readInput(in_file_name):
    """Read the input from files.
    Read maternal alleles, paternal alleles, and samples per position.
    
    Arguments:
    file name -- string
    returns lists of parsed data    
    
    """
    in_file = open(in_file_name, 'r')
    positions = []
    samples = []
    M = []; P = [];
    MC = []; PC = [];
    while True:
        line = in_file.readline()
        if not line: break
        if line[0] == '#': continue  #skip comment
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
            binaryWeights = map(float, line[1].split())
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


def computeEval(reference, prediction, sensitivity, num_patt):
    """Compute evaluation of referenced anotation versus predicted one.
    
    Prediction for one continuous region R of an inheritance pattern is considered
    correct if there is a continuous part with length at least '1/sensitivity'*len(R)
    that is predicted correctly.
    
    """
    num_real = 0
    num_called = 0
    ok_called = []
    not_called = []
    
    ind = 0
    while ind < len(reference):
        ctype = reference[ind]
        length = 0
        begin_pos = ind
        while ind < len(reference) and reference[ind] == ctype:
            length += 1
            ind += 1
        end_pos = ind
        
        max_ok_length = 0
        tmp_length = 0
        stats = [0 for x in range(num_patt)]
        for i in range(begin_pos, end_pos):
            stats[prediction[i]] += 1
            if prediction[i] == ctype:
                tmp_length += 1
            else:
                #print i, prediction[i],  ctype, type(prediction[i]), type(ctype)
                max_ok_length = max(max_ok_length, tmp_length)
                tmp_length = 0
        max_ok_length = max(max_ok_length, tmp_length)
        
        num_real += 1
        if max_ok_length >= length/float(sensitivity):
            num_called += 1
            ok_called.append((length, ctype))
        else:
            max_predic = -1
            max_type = -1
            for i in range(len(stats)):
                if stats[i] > max_predic:
                    max_predic = stats[i]
                    max_type = i
            not_called.append((length, ctype, max_type))
        #else: print  max_ok_length, length, reference[ind], prediction[ind], ctype
    return num_called, num_real, ok_called, not_called

def categorize(data, by_length, by_type):
    lengths = [100, 500, 1000, 5000, 10000, 20000, 50000, 100000]
    for item in data:
        length = item[0]
        ptype = item[1]
        gen_len = sorted([(abs(length-l), l) for l in lengths])[0][1]
        #print length, gen_len
        
        #by length
        if not gen_len in by_length: by_length[gen_len] = 0
        if ptype != 3: by_length[gen_len] += 1 #don't count normal IP into CNV length stats
        #by type
        if not ptype in by_type: by_type[ptype] = 0
        by_type[ptype] += 1

def computeStats(ok, wrong, pref, num_patt):
    lengths = [100, 500, 1000, 5000, 10000, 20000, 50000, 100000]
    ok_by_length = dict([(i,0) for i in lengths])
    ok_by_type = dict([(i,0) for i in range(num_patt)])
    wr_by_length = dict([(i,0) for i in lengths])
    wr_by_type = dict([(i,0) for i in range(num_patt)])
    categorize(ok, ok_by_length, ok_by_type)
    categorize(wrong, wr_by_length, wr_by_type)
    for l in lengths:
        o = ok_by_length[l]
        w = wr_by_length[l]
        ratio = 0
        if o+w > 0: ratio = o/float(o+w) *100
        print pref, l, ": ", o, w, ratio, '%'
    for t in range(num_patt):
        o = ok_by_type[t]
        w = wr_by_type[t]
        ratio = 0
        if o+w > 0: ratio = o/float(o+w) *100
        print pref, t, ": ", o, w, ratio, '%'

def test(fcnv, snp_positions, samples, M, P, MSC, PSC, mixture, ground_truth, file_name_prefix):
    vp, v_state_path = fcnv.viterbiPath(samples, M, P, MSC, PSC, mixture) 
    #vp = fcnv.maxPosteriorDecoding(samples, M, P, MSC, PSC, mixture)
    
    date = datetime.now().strftime('%m-%d-%H-%M')
    #fout = file(file_name_prefix + ".prediction" + date + ".txt", 'w')
    annot_out = file(file_name_prefix + ".annotation" + date + ".txt", 'w')
    
    
    print fcnv.inheritance_patterns
    num_patt = fcnv.getNumIP()

    viterbi_correct = 0
    state_correct_vp = [0 for x in range(num_patt)]

    for i in xrange(len(vp)):
        print >>annot_out, snp_positions[i], ground_truth[i], vp[i] #, v_state_path[i]

        viterbi_correct += int(ground_truth[i] == vp[i])
        state_correct_vp[ground_truth[i]] += int(ground_truth[i] == vp[i])

    print "viterbi: ", state_correct_vp

    print 'Viterbi  : ',(viterbi_correct*100.)/len(vp), '% OK'

    if len(vp) != len(ground_truth):
        print "UNEXPECTED ERROR: different prediction lengths:", len(vp), len(ground_truth)
        return
    
    for i in [2]: #, 2, 4, 8, 16]:
        #precision and recall
        print "sensitivity: 1 /", i
        
        #recall Viterbi
        called_v, real, ok_called, not_called = computeEval(ground_truth, vp, i, num_patt)
        print "viterbi recall   : ", called_v, '/',  real, " => ", 
        print "%0.3f" % (called_v*100./real), "%"
        print ok_called
        print not_called
        computeStats(ok_called, not_called, "VR", num_patt)
        
        #precision Viterbi 
        correct_v, claimed_v, ok_prediction, wr_prediction = computeEval(vp, ground_truth, i, num_patt)
        print "viterbi precision: ", correct_v , '/',  claimed_v, " => ", 
        print "%0.3f" % (correct_v*100./claimed_v), "%"
        print ok_prediction
        print wr_prediction
        computeStats(ok_prediction, wr_prediction, "VP", num_patt)
        
        print "%0.3f" % (called_v*100./real), "\%"
        print "%0.3f" % (correct_v*100./claimed_v), "\%"
        
        
def main():
    #parse command line arguments
    parser = argparse.ArgumentParser(description='Performs fetal CNV analysis from maternal plasma and phased parental data.')
    parser.add_argument('input', type=str, nargs=1, help='path to input file with allele counts in plasma and parental haplotypes')
    parser.add_argument('target', type=str, nargs=1, help='path to file with background truth - "target"')
    parser.add_argument('plasma', type=str, nargs=1, help='path to file with plasma sequencing DOC for all chromosomal positions')
    parser.add_argument('ref', type=str, nargs=1, help='path to file with reference plasma sequencing DOC for all chromosomal positions')
    parser.add_argument('seq', type=str, nargs=1, help='path to ref. genomic sequence in fasta format')
    parser.add_argument('param', type=str, nargs=1, help='path to file with method parameters')
    parser.add_argument('--ff', type=float, help='fetal mixture ratio', default=-1.)
    parser.add_argument('--useCvrg', help='use coverage flag', action="store_true")
    parser.add_argument('--trainGrad', type=str, help='train by maxll gradient and output new params to given file', default="")
    parser.add_argument('--getObsCounts', help='get observed allele counts', action="store_true")
    args = parser.parse_args()
    
    in_file_name = args.input[0]
    target_file_name = args.target[0]
    plasma_doc_file = open(args.plasma[0], "r")
    ref_doc_file = open(args.ref[0], "r")
    seq_file = open(args.seq[0], "r")
    param_file = open(args.param[0], "r")
    if args.ff > 0: mix = args.ff
    runGradTraining = False
    if args.trainGrad != "":
         runGradTraining = True
         res_param_file_name = args.trainGrad
    
    #print input info
    print "------------------------------------------"
    print "Running fCNV, input parameters:"
    print "input:", in_file_name
    print "target:", target_file_name
    print "plasma:", plasma_doc_file
    print "refDOC:", ref_doc_file
    print "seq:", seq_file
    print "param:", param_file
    print "--ff:", args.ff
    print "--useCvrg:", args.useCvrg
    print "--trainGrad:", args.trainGrad
    print "------------------------------------------"
    os.system("hostname")
    
    #read the pre-processed input
    snp_positions, samples, M, P, MSC, PSC = readInput(in_file_name)
    
    #fetch the method parameters from file
    crfParams = readParams(param_file)
    print "============ CRF PARAMETERS =============="
    for p in sorted(crfParams):
        print p, "=", crfParams[p]
    print "=========================================="
    
    #get genomic positions on the last lines of the pileup files to estimate the length of the chromosome
    if args.useCvrg:
        with open(args.plasma[0], 'rb') as fh:
            fh.seek(-256, 2)
            last_pos_plasma = int(fh.readlines()[-1].decode().split(' ')[0])
            fh.close()
#        with open(args.ref[0], 'rb') as fh:
#            fh.seek(-256, 2)
#            last_pos_ref = int(fh.readlines()[-1].decode().split(' ')[0])
#            fh.close()
        chr_length = last_pos_plasma + 4742
        
        gc_sum = [0] * chr_length
        prefix_sum_plasma = [0] * chr_length
        prefix_count_plasma = [0] * chr_length
        prefix_sum_ref = [0] * chr_length
        prefix_count_ref = [0] * chr_length
    
        #get GC content prefix sums from the reference
        gen_pos = 0
        keep_reading = True
        while keep_reading:
            line = seq_file.readline().strip().upper()
            if len(line) == 0: break
            if line[0] == '>': continue
            for i in range(len(line)):
                gc_sum[gen_pos] = gc_sum[max(gen_pos - 1, 0)]
                if line[i] in 'GC': gc_sum[gen_pos] += 1
                gen_pos += 1
                if gen_pos >= chr_length:
                    keep_reading = False
                    break
        seq_file.close()
        
        last = 0
        while True: 
            line = plasma_doc_file.readline()
            if not line: break
            row = map(int, line.split(' '))
            if row[0] >= chr_length: break
            prefix_sum_plasma[row[0]] = prefix_sum_plasma[last] + row[1]
            prefix_count_plasma[row[0]] = prefix_count_plasma[last] + 1
            last = row[0]
        plasma_doc_file.close()
        
        last = 0
        while True: 
            line = ref_doc_file.readline()
            if not line: break
            row = map(int, line.split(' '))
            if row[0] >= chr_length: break
            prefix_sum_ref[row[0]] = prefix_sum_ref[last] + row[1]
            prefix_count_ref[row[0]] = prefix_count_ref[last] + 1
            last = row[0]
        ref_doc_file.close()
    #ENDIF 
    
    #snp_positions = []
    ground_truth = []
    target_file = open(target_file_name, 'r')
    while True: 
        line = target_file.readline()
        if not line: break
        line = line.rstrip("\n").split("\t")
        #snp_positions.append(int(line[0]))
        ground_truth.append(int(line[-1]))
    target_file.close()
    
    fcnv = fcnvCRF.FCNV(crfParams, None, None, False)
    mix, mix_median, ct = fcnv.estimateMixture(samples, M, P)
    print "Est. Mixture: ", mix, mix_median, '(', ct ,')',
    #mix = 0.13 #proportion of fetal genome in plasma
    if args.ff > 0: mix = args.ff
    print "used:", mix
    
    cnv_prior = None
    if args.useCvrg:
        pass
#        cvrg = cvrgHMM.coverageFCNV(snp_positions, prefix_sum_plasma, prefix_count_plasma, prefix_sum_ref, prefix_count_ref, gc_sum)
#        cvrg_posterior = cvrg.posteriorDecoding(mix)
#        del prefix_sum_plasma, prefix_count_plasma, prefix_sum_ref, prefix_count_ref, gc_sum
#        
#        cnv_prior = [ [0., 0., 0.] for x in range(len(snp_positions)) ]
#        for pos in range(len(cvrg_posterior)):
#            for cp_num_posterior in cvrg_posterior[pos]:
#                cnv_prior[pos][cp_num_posterior[1]] = cp_num_posterior[0]
#        del cvrg, cvrg_posterior
    
    del fcnv
    fcnv = fcnvCRF.FCNV(crfParams, snp_positions, cnv_prior, args.useCvrg)
    
    ground_truth = []
    target_file = open(target_file_name, 'r')
    for line in target_file.readlines():
        line = line.rstrip("\n").split("\t")
        ground_truth.append(int(line[-1]))
    target_file.close()
    
    parameterStats = dict()
    
    
    #get get observed allele counts
    if args.getObsCounts:
        parameterStats = fcnv.computeCountsTable(ground_truth, samples, M, P, MSC, PSC, mix, parameterStats)
        varFileName = target_file_name.split('/')[-1].split('.')[0].replace(':', '-')+'.observedCounts.txt'
        varFile = open(varFileName, 'w')
        print >>varFile, len(parameterStats)
        for expP in parameterStats.keys():
            print >>varFile, expP
            print >>varFile, " ".join(map( lambda x: '/'.join(map(str, x)), parameterStats[expP]) )
        varFile.close()
        return 0
    
    #run gradient training
    if runGradTraining:
        #run the training iterations
        for iterNum in range(1):
        #    print "iterNum: ", iterNum
            #ll, params = fcnv.computeLLandGradient(ground_truth, samples, M, P, MSC, PSC, mix) 
            pregts, preps, preloss, params, postgts, postps, postloss = fcnv.computeLLandMaxMarginUpdate(ground_truth, samples, M, P, MSC, PSC, mix, float("inf"), compute_postloss=True)
            print preloss, params
            print "{0} !>= {1}".format(pregts - preps, preloss)
            print "{0} >= {1}".format(postgts - postps, postloss)
            print "------------------------------------------------------------------"
        
        #save the trained parameters to the file
        res_param_file = open(res_param_file_name, "w")    
        for p in sorted(params):
            if isinstance(params[p], list):
                print >>res_param_file, p, "=", " ".join(map(str, params[p]))
            else:
                print >>res_param_file, p, "=", params[p]
        res_param_file.close()
        
        return 0
    
    #res_file = open(out_file_name, 'w')
    file_name_prefix = target_file_name.split('/')[-1].split('.')[0].replace(':', '-')
    print "------------------ w/o TRAINING -------------------"
    test(fcnv, snp_positions, samples, M, P, MSC, PSC, mix, ground_truth, file_name_prefix)
    #test(fcnv, snp_positions[:1000], samples[:1000], M[:1000], P[:1000], MSC[:1000], PSC[:1000], mix, ground_truth[:1000], file_name_prefix)
    
    
if __name__ == "__main__":
    #import doctest
    #doctest.testmod()
    
    main()
    



