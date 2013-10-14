#!/usr/bin/pypy
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import argparse

k=1000
binSizes = [0, 100, 300, 500, 1000, 5000, 10000] 
numBins = 7

truePos= {'MDup': [0]*numBins, 'MDel': [0]*numBins,'PDup': [0]*numBins,'PDel': [0]*numBins, 'sum':[0]*numBins}
falsePos={'MDup': [0]*numBins, 'MDel': [0]*numBins,'PDup': [0]*numBins,'PDel': [0]*numBins, 'sum':[0]*numBins} 
total= {'MDup': [0]*numBins, 'MDel': [0]*numBins,'PDup': [0]*numBins,'PDel': [0]*numBins, 'sum':[0]*numBins}

binSizes= map(lambda x: x*k, binSizes)

def findIndex(sz):
    difs=map(lambda x: abs(sz-x), binSizes)
    return difs.index(min(difs))

def main():
    parser = argparse.ArgumentParser(description='This script filters the given reads of a certain region of plasma so that it simulates a deletion in this region. It requires the read file, snips file, the rate of fetus DNA in the plasma and the target haplotype of fetus for deletion')
    parser.add_argument('summaryFile', type=str, nargs=1, help='the address to the summary file')
    args = parser.parse_args()

    summary_file=open(args.summaryFile[0], "r")
   
    nameMap={}
    
    for line in summary_file:
        if line.find(".txt")!=-1:
            rawHeader=line
        if line.find('Real State:')!=-1:
            parts=line.split(' ')
            if (rawHeader.find('cvrg')==-1 and int(parts[2])==2) or (rawHeader.find('cvrg')!=-1 and int(parts[2])==0):
                nameMap[rawHeader]='PDel'
            if (rawHeader.find('cvrg')==-1 and int(parts[2])==0):
                nameMap[rawHeader]='MDel'
            if (rawHeader.find('cvrg')==-1 and int(parts[2])==6) or (rawHeader.find('cvrg')!=-1 and int(parts[2])==2):
                nameMap[rawHeader]='MDup'
            if (rawHeader.find('cvrg')==-1 and int(parts[2])==4):
                nameMap[rawHeader]='PDup'
    



    summary_file.seek(0)
    for line in summary_file:
        
        #if a new result has started
        if line.find(".txt")!=-1:
            rawHeader=line
            header= line.split('-')
            if (header[4].startswith('delete')):
                simSZ=int(header[3])-int(header[2])
            else:
                simSZ=int(header[4])-int(header[3])
        if line.find('predicted')!=-1:
            parts=line.split(' ')
            reg= map(int, parts[0].split('-'))
            sz=reg[1]-reg[0]
            if int(parts[2])!=int(parts[5]) and ((rawHeader.find('cvrg')==-1 and (int(parts[5])==3)) or (rawHeader.find('cvrg')!=-1 and int(parts[5])==1)):
                falsePos[nameMap[rawHeader]][findIndex(sz)]+=1
                falsePos['sum'][findIndex(sz)]+=1
        if (line.startswith('Recall')):
            truePos[nameMap[rawHeader]][findIndex(simSZ)]+= int(line.split('\t')[1])
            truePos['sum'][findIndex(simSZ)]+= int(line.split('\t')[1])
        if (line.startswith('Recall')):
            total[nameMap[rawHeader]][findIndex(simSZ)]+=1
            total['sum'][findIndex(simSZ)]+=1
#        if (line.startswith('Recall')):
#            if int(line.split('\t')[1])==0 and simSZ==1000*k:
#                print rawHeader


    print '*',
    print "\t",
    print 'Size :',
    print "\t",
    print "\t",
    for (i,num) in enumerate(binSizes):
            print num/k,
            print "k\t",

    print 'total',
    for key in sorted(total.keys()):
        print ''
        print key,
        print '\tRecall :',
        print "\t",
        for (i,num) in enumerate(binSizes):
            if total[key][i]==0:
                print 1,
                print "\t",
            else: 
                print '%0.4f' % (float(truePos[key][i])/total[key][i]),
                print "\t",
        #        print "Recall for ",
        #        print num,
        #        print ":: ",
        if sum(total[key])==0:
            print 1,
            print "\t",
        else:
            print '%0.4f' % (float(sum(truePos[key]))/sum(total[key])),
            print "\t",

        print ''
        print key,
        print '\tPrecision :',
        print "\t",
        for (i,num) in enumerate(binSizes):
            #print float(truePos[i])/(truePos[i]+falsePos[i]),
            print str(truePos[key][i])+'/'+str(truePos[key][i]+falsePos[key][i]),
            print "\t",
        print str(sum(truePos[key]))+'/'+str(sum(truePos[key])+sum(falsePos[key])),
        print "\t",
        print ''
    #        print "Precision for ",
    #        print num,
    #        print ":: ",



 

if __name__ == '__main__':
    main()

