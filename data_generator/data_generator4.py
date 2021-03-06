import random
import math
import sys

def allelesDistributionMeans(maternal, fetal, mix):
        """        
        Arguments:
        maternal -- maternal alleles
        paternal -- paternal alleles
        mix -- fetal admixture ratio
        returns list(floats) -- list of nucleotides probabilities
        """
        adjusted_fetal_admix = mix/2. * len(fetal)
        adjusted_maternal_admix = (1.-mix) / 2. * len(maternal)
        cmix = adjusted_fetal_admix / (adjusted_maternal_admix + adjusted_fetal_admix)
        dist = []
        for nuc in nucleotides:
            p = maternal.count(nuc) / float(len(maternal)) * (1.-mix)
            p += fetal.count(nuc) / float(len(fetal)) * adjusted_fetal_admix
            dist.append(p)

        return dist
        
#parse command line arguments => suffix of the output files
suffix = ""
if len(sys.argv) == 2: suffix = sys.argv[1]
    

nucleotides = ['A', 'C', 'G', 'T']

length = 10**5 #how many "SNP" positions to generate
mixture = 0.10 #proportion of fetal genome in plasma
mean_coverage = 80 #average coverage of the plasma \mu
#dev_coverage = math.sqrt(mean_coverage)  #deviation \sigma
haptypeblock_mu = 100
haptypeblock_dev = 20

#how many SNP an CNV event will span - normal distribution
#cnv_length_mu = 1000
#cnv_length_dev = 200
cnv_length_mus = [30, 60] #[100, 500, 1000, 5000]#, 10000, 20000, 50000]
cnv_length_devs = [math.sqrt(x) for x in cnv_length_mus] #deviation is sqrt of the length

file_paternal = "paternal_alleles" + suffix + ".txt"
file_maternal = "maternal_alleles" + suffix + ".txt"
file_fetal = "fetal_alleles" + suffix + ".txt"
file_samples = "plasma_samples" + suffix + ".txt"

fpaternal = open(file_paternal, 'w')
fmaternal = open(file_maternal, 'w')
ffetal = open(file_fetal, 'w')
fsamples = open(file_samples, 'w')

#generate the haplotypes
M = []
P = []
for i in range(length):
    while True:
        alleles_set = []
        alleles_set.append(random.choice(nucleotides))
        alleles_set.append(random.choice(nucleotides))
        if alleles_set[0] != alleles_set[1]: break
    mnew = []
    pnew = []
    mnew.append(random.choice(alleles_set))
    mnew.append(random.choice(alleles_set)) 
    pnew.append(random.choice(alleles_set))
    pnew.append(random.choice(alleles_set))
    
    #print mnew, pnew, "=========="
    M.append(mnew)
    P.append(pnew)
    print >>fmaternal, " ".join(mnew)
    print >>fpaternal, " ".join(pnew)
fmaternal.close()
fpaternal.close()

#generate inheritance patterns
inheritance_patterns = []
for mats in range(3):
    for pats in range(3):
        if mats+pats in [0,4]: continue
        inheritance_patterns.append((mats, pats)) 
num_patterns = len(inheritance_patterns)
        
#generate the map of inheritance patterns in the generated sequence
pattern_map = [3 for i in xrange(length)]
#for patt in range(len(inheritance_patterns)):
#    if patt == 3: continue
#    for size_id in range(len(cnv_length_mus)):
#        #choose how many SNP the event will span
#        size = int(random.gauss(cnv_length_mus[size_id], cnv_length_devs[size_id]))
#        
#        #pick a position and check if the space is available
#        start_pos = -1
#        while True:
#            start_pos = random.randint(0, length-1)
#            if start_pos + size > length: continue
#            for pos in range(start_pos, start_pos + size):
#                if not pattern_map[pos] == 3: continue
#            break;
#        
#        #write the CNV to the pattern_map
#        for pos in range(start_pos, start_pos + size):
#            pattern_map[pos] = patt
            
#initialize hidden haplotype variables
#first two fields indicate haplotypes for 2 possible maternal copies
# the other two the same for paternal
htypes = [0, 0, 0, 0]
hblock_lengths = [0, 0, 0, 0]
#generate fetal alleles
F = []
for pos in range(length):
    #inherited haplotype blocks
    for i in range(len(hblock_lengths)):
        hblock_lengths[i] -= 1
        if hblock_lengths[i] <= 0:
          hblock_lengths[i] = int(round(random.gauss(haptypeblock_mu, haptypeblock_dev)))
          htypes[i] = random.choice([0,1])
    
    #get the type of the event (inheritence pattern)
    X = pattern_map[pos]
    
    #generate alleles for this position
    fnew = []
    for k in range(inheritance_patterns[X][0]):
        hap = htypes[0+k]
        noise = int(random.random() >= 0.9)
        hap = (hap + noise) % 2
        fnew.append( M[pos][hap] )
    for k in range(inheritance_patterns[X][1]):
        hap = htypes[2+k]
        noise = int(random.random() >= 0.9)
        hap = (hap + noise) % 2
        fnew.append( P[pos][hap] )
        
    F.append(fnew)
    print >>ffetal, " ".join(fnew), X
ffetal.close()

#generate "reads"
for i in range(length):
    ps = allelesDistributionMeans(M[i], F[i], mixture)
    samples = [ max(0, int(round(random.gauss(mean_coverage*ps[x], max(0.7, math.sqrt(mean_coverage*ps[x])))))) for x in range(4) ]
    print >>fsamples, " ".join(map(str, samples))
fsamples.close()   




