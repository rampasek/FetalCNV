#!/usr/bin/python2

import argparse
import subprocess

prev_lines = ['', '']
def getlineFromFile(in_files, out_files, idn):
    global prev_lines
    
    if prev_lines[idn] != '':
        print >>out_files[idn], prev_lines[idn],
    
    line = in_files[idn].readline()
    while len(line) > 0 and line[0] == '#':
        #ignore the header
        print >>out_files[idn], line,
        line = in_files[idn].readline()
    prev_lines[idn] = line
    
    return line
    

def main():
    global prev_lines
    #parse ARGs
    parser = argparse.ArgumentParser(description='Filter SNP positions by calling quality and min. coverage. Awaits filenames for M, and P .vcf files.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to .vcf files with M, and P SNPs')
    args = parser.parse_args()
    
    if len(args.filenames) != 2: die("Unexpected number of arguments passed! Expecting 2 filenames.")
    
    #treat these as CONSTANTS!
    M = 0; P = 1; 
    ALL = [M, P]
    
    #list of input files
    in_files = [open(args.filenames[i], "r" ) for i in ALL]
    
    #list of output files
    out_files = [open(args.filenames[i][:-3]+"ftr.vcf", "w") for i in ALL]
    
    #positions ignored from the union of M and P
    ignored_pos = 0
    
    #read SNPs from M, P, F vcf files
    snps = [[] for i in ALL]
    for i in ALL:
        snps[i] = getlineFromFile(in_files, out_files, i).split('\t')
    
    #output genotypes for all positions in UNION of M and P SNP positions
    while len(snps[M])>2 or len(snps[P])>2: #while there is a SNP positions in M or P
        #get the position of SNP that occure first
        for i in ALL:
            if snps[i][0] == '': #if an input files is already at EOF
                snps[i][0] = 'chrZZ'
                snps[i].append(1e15)
            else: #convert to int
                snps[i][1] = int(snps[i][1])
        #chromosome
        min_chr = min(snps[M][0], snps[P][0])
        #position
        min_pos = 1e15
        for i in ALL:
            if min_chr == snps[i][0] and snps[i][1] < min_pos:
                min_pos = snps[i][1]
        
        callQ = [0. for i in ALL]
        for i in ALL:
            #if there is a SNP in the data at this position
            if min_chr == snps[i][0] and min_pos == snps[i][1]:
                info = snps[i]
                if len(info) <= 2: continue
                #parse out quality info
                callQ[i] = float(info[5])
        qualityOK = bool(callQ[M] >= 50 or callQ[P] >= 50)
        
        coverage = [0 for i in ALL]
        for i in ALL:
            cmd = 'samtools mpileup -r %(chr)s:%(pos)d-%(pos)d __%(gnm)s.part.bam' % {'chr':min_chr, 'pos':min_pos, 'gnm':'MP'[i]}
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = process.communicate()
            if process.returncode == 0:
                fields = out.split('\t')
                if len(fields) < 4:
                    #there is no coverage in the BAM file for this position
                    coverage[i] = 0
                    #print "!!! :", out, "|", err, "|", cmd
                else:
                    coverage[i] = int(fields[3])
                #print out, '>>coverage>>', coverage[i]
            else:
                print err
        coverageOK = bool(coverage[M] >= 15 and coverage[P] >= 15)
        
        MOK = bool(callQ[M] >= 50 or coverage[M] >= 15)
        POK = bool(callQ[P] >= 50 or coverage[P] >= 15)
        if not (MOK and POK): #ignore positions that are not good enough
            ignored_pos += 1
            for i in ALL:
                #if there is a SNP in the data at this position, skip it
                if min_chr == snps[i][0] and min_pos == snps[i][1]:
                    #print min_pos, callQ[M], callQ[P], coverage[M], coverage[P],  prev_lines[i],
                    prev_lines[i] = ''
                    
            
        
        #read input: next SNP
        for i in ALL:
            if min_chr >= snps[i][0] and min_pos >= snps[i][1]:
                snps[i] = getlineFromFile(in_files, out_files, i).split('\t')
                
   print "Low quality positions ignored in the region:", ignored_pos  
    
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    
    main()

