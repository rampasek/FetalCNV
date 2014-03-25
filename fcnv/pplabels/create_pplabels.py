#!/usr/bin/pypy

import argparse
import sys
import random
import copy
import itertools
from datetime import datetime

class HMMState(object):
    """Class to represent a particular phased inheritance pattern."""
    
    def __init__(self, ip, pp):
        """Initialize new HMMState object
        Arguments:
        ip -- tuple(int, int) -- general inheritance pattern
        pp -- tuple(tuple(ints), tuple(ints)) -- phased inheritance pattern
        """
        super(HMMState, self).__init__()
        
        self.inheritance_pattern = ip
        self.phased_pattern = pp
        
    def __str__(self):
        return "IP: {0}, PP: {1}".format(self.inheritance_pattern, self.phased_pattern)
        
    def __eq__(self, other):
        return self.inheritance_pattern == other.inheritance_pattern and \
            self.phased_pattern == other.phased_pattern

def main():
    #parse command line arguments
    parser = argparse.ArgumentParser(description='Create phased pattern target from inheritance pattern target file using child\'s haplotypes annotation')
    parser.add_argument('chHapsAnnot', type=str, nargs=1, help='path to file with child\'s haps\' origin annotation')
    parser.add_argument('target', type=str, nargs=1, help='path to file with unphased IP background truth - "target"')
    parser.add_argument('pplabels', type=str, nargs=1, help='name of output file with phased pattern target labels')
    args = parser.parse_args()
    
    ### read the unphased IP target file
    target_file = open(args.target[0], "r")
    target = []
    for line in target_file.readlines():
        line = line.rstrip("\n").split("\t")
        line[0] = int(line[0])
        line[-1] = int(line[-1])
        target.append(line)
    target_file.close()
    #for x in target: print x
    
    
    ### child's snp positions and their haloptype annotation
    # -- particular maternal/paternal haplotype that got inherited
    channot_file = open(args.chHapsAnnot[0], "r")
    snp_positions = []
    chAnnot = []
    for line in channot_file.readlines():
        line = map(int, line.rstrip('\n').split())
        snp_positions.append(line[0])
        chAnnot.append(line[1])
    channot_file.close()
    
    #check input 
    #for i in range(len(snp_positions)):
    #    print snp_positions[i], chAnnot[i]
    
    ### generate the same patterns as used in the chHapsAnnot input
    ch_haplo_patterns = []
    ip = (1, 1)
    for mPhased in itertools.combinations_with_replacement([0,1], ip[0]):
        for pPhased in itertools.combinations_with_replacement([0,1], ip[1]):
            pattern = (tuple(['M'] + list(mPhased)), tuple(['P'] + list(pPhased)))
            ch_haplo_patterns.append(pattern)
    for pPhased in itertools.combinations_with_replacement([0,1], ip[0]):
        for mPhased in itertools.combinations_with_replacement([0,1], ip[1]):
            pattern = (tuple(['P'] + list(pPhased)), tuple(['M'] + list(mPhased)))
            ch_haplo_patterns.append(pattern)
    #for i, patt in enumerate(ch_haplo_patterns): print i, patt
    
    ### generate the phased_patterns as used in CRF/HMM
    #generate inheritance patterns
    inheritance_patterns = []
    for mats in range(3):
        for pats in range(3):
            if mats+pats in [0,4]: continue
            inheritance_patterns.append((mats, pats))
 
    #generate phased variants of the inheritance patterns
    states = []
    phased_patterns = []
    for ip in inheritance_patterns:
        num_phased_states = 0
        for mPhased in itertools.combinations_with_replacement([0,1], ip[0]):
            for pPhased in itertools.combinations_with_replacement([0,1], ip[1]):
                phased_pattern = (tuple(mPhased), tuple(pPhased))
                phased_patterns.append(phased_pattern)
                new_state = HMMState(ip, phased_pattern)
                states.append(new_state)
    #for i, state in enumerate(states): print i, state
    
    
    ### parse CNV type of the data-case from the target file name
    try:
        fname = args.target[0].split('/')[-1]
        parts = fname.split('-')
        
        duplication = False
        deletion = False
        hap_origin = -1
        if fname.find('duplicate') != -1:
            duplication = True
            hap_origin = ['a', 'b'].index(parts[1].lower())
        if fname.find('delete') != -1:
            deletion = True
            hap_origin = ['a', 'b'].index(parts[0].lower())
        
        maternal_origin = False
        paternal_origin = False
        if parts[0].lower().find('m') != -1:
            maternal_origin = True
        elif parts[0].lower().find('p') != -1:
            paternal_origin = True
    except Exception as e:
        print e, fname
    
    print 'dup='+str(duplication), 'del='+str(deletion), 'hap='+str(hap_origin), 'matO='+str(maternal_origin), 'patO='+str(paternal_origin)
 
#    0 (('M', 0), ('P', 0))  --> 7 phased_pattern
#    1 (('M', 0), ('P', 1))  --> 8 phased_pattern
#    2 (('M', 1), ('P', 0))  --> 9 phased_pattern
#    3 (('M', 1), ('P', 1))  --> 10 phased_pattern
#    4 (('P', 0), ('M', 0))  --> 7 phased_pattern
#    5 (('P', 0), ('M', 1))  --> 9 phased_pattern
#    6 (('P', 1), ('M', 0))  --> 8 phased_pattern
#    7 (('P', 1), ('M', 1))  --> 10 phased_pattern
    chAnnot2ppIndex = [7, 8, 9, 10, 7, 9, 8, 10]

#        0       1       2       3       4       5       6
#    [(0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1)]
#    0  IP: (0, 1), PP: ((),     (0,)  )
#    1  IP: (0, 1), PP: ((),     (1,)  )
#    2  IP: (0, 2), PP: ((),     (0, 0))
#    3  IP: (0, 2), PP: ((),     (0, 1))
#    4  IP: (0, 2), PP: ((),     (1, 1))
#    5  IP: (1, 0), PP: ((0,),   ()    )
#    6  IP: (1, 0), PP: ((1,),   ()    )
#    7  IP: (1, 1), PP: ((0,),   (0,)  )
#    8  IP: (1, 1), PP: ((0,),   (1,)  )
#    9  IP: (1, 1), PP: ((1,),   (0,)  )
#    10 IP: (1, 1), PP: ((1,),   (1,)  )
#    11 IP: (1, 2), PP: ((0,),   (0, 0))
#    12 IP: (1, 2), PP: ((0,),   (0, 1))
#    13 IP: (1, 2), PP: ((0,),   (1, 1))
#    14 IP: (1, 2), PP: ((1,),   (0, 0))
#    15 IP: (1, 2), PP: ((1,),   (0, 1))
#    16 IP: (1, 2), PP: ((1,),   (1, 1))
#    17 IP: (2, 0), PP: ((0, 0), ()    )
#    18 IP: (2, 0), PP: ((0, 1), ()    )
#    19 IP: (2, 0), PP: ((1, 1), ()    )
#    20 IP: (2, 1), PP: ((0, 0), (0,)  )
#    21 IP: (2, 1), PP: ((0, 0), (1,)  )
#    22 IP: (2, 1), PP: ((0, 1), (0,)  )
#    23 IP: (2, 1), PP: ((0, 1), (1,)  )
#    24 IP: (2, 1), PP: ((1, 1), (0,)  )
#    25 IP: (2, 1), PP: ((1, 1), (1,)  )

    
    ### create phased pattern target labeling
    out_file_name = args.pplabels[0]
    if out_file_name == 'auto':
        out_file_name = args.target[0].replace('target', 'pplabels')
    #print out_file_name
    
    out_file = open(out_file_name, "w")
    offset = snp_positions.index(target[0][0])
    #print "offset", offset
    for i in range(len(target)):
        while target[i][0] > snp_positions[offset+i]:
            offset +=1
        if target[i][0] != snp_positions[offset+i]:
            print "ERROR: SNP position doesn't match", target[i][0], snp_positions[offset+i];
            return -1
        
        
        chInheritance = chAnnot2ppIndex[chAnnot[offset+i]]
        ip = target[i][-1]
        pptarget = -1
        
        #Normal
        if   ip == 3: 
            pptarget = chInheritance
        
        #Maternal Deletion
        elif ip == 0:
            if not deletion:
                print "ERROR: MDel", target[i]
                return -1
            #maternal part of what the child inherited was deleted
            if chInheritance in [7, 9]:
                pptarget = 0
            elif chInheritance in [8, 10]:
                pptarget = 1
        
        #Paternal Deletion
        elif ip == 2:
            if not deletion:
                print "ERROR: PDel", target[i]
                return -1
            #paternal part of what the child inherited was deleted
            if chInheritance in [7, 8]:
                pptarget = 5
            elif chInheritance in [9, 10]:
                pptarget = 6
            
        #Maternal Duplication
        elif ip == 6:
            if not (maternal_origin and duplication and hap_origin != -1):
                print "ERROR: MDup", target[i], hap_origin
                return -1
                
            patt = list(phased_patterns[chInheritance])
            patt = map(list, patt)
            patt[0].append(hap_origin)
            patt[0].sort()
            patt = map(tuple, patt)
            patt = tuple(patt)
            
            pptarget = phased_patterns.index(patt)
            if not ( 20 <= pptarget <= 25):
                print "ERROR: MDup: pattern index", pptarget
                return -1
        
        #Paternal Duplication
        elif ip == 4:
            if not (paternal_origin and duplication and hap_origin != -1):
                print "ERROR: PDup", target[i], hap_origin
                return -1
                
            patt = list(phased_patterns[chInheritance])
            patt = map(list, patt)
            patt[1].append(hap_origin)
            patt[1].sort()
            patt = map(tuple, patt)
            patt = tuple(patt)
            
            pptarget = phased_patterns.index(patt)
            if not ( 11 <= pptarget <= 16):
                print "ERROR: PDup: pattern index", pptarget
                return -1
        
        else:
            print "ERROR: unexpected IP id in the target file:", ip
            return -1
        
        if pptarget == -1:
            print "ERROR: new phased pattern target not determined for", target[i]
            return -1
        
        target[i][-1] = pptarget
        print  >>out_file, '\t'.join(map(str, target[i])) #, pptarget
    out_file.close()
    print 'DONE', out_file_name
    
    
if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    
    main()

