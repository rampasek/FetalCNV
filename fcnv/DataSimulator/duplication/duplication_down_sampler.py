#!/usr/bin/python2

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
import random

def main():
    parser = argparse.ArgumentParser(description='Filter SNP positions by call quality and min. coverage. Awaits filenames for M, P .vcf files, and M, P .sam files.')
    parser.add_argument('filenames', type=str, nargs='+', help='paths to .vcf files with M, P SNPs and to corresponding .sam files')
    args = parser.parse_args()

    readsFiltered_file= open(args.filenames[0], "r")

    targetCoverage= float(args.filenames[1]) * float(args.filenames[2]) / 2.0

    region= args.filenames[5].strip().split(':')[1].strip().split('-')
    currentCoverage= float(args.filenames[3]) * int(args.filenames[4]) / (int(region[1])-int(region[0]))

    outputRate= targetCoverage/currentCoverage

    #print targetCoverage
    #print currentCoverage
    #print outputRate
    #exit()

    for line in readsFiltered_file:
        if random.random() < outputRate:             
                print line,

if __name__ == '__main__':
    main()

