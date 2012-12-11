import random
import math

nucleotides = ['A', 'C', 'G', 'T']

length = 10**5 #how many "SNP" positions to generate
mixture = 0.1 #proportion of fetal genome in plasma
mean_coverage = 75 #average coverage of the plasma \mu
dev_coverage = 5   #deviation \sigma

#how many SNP an CNV event will span - uniform distrib
#cnv_length_min = 500
#cnv_length_max = 1500
cnv_length_mu = 300
cnv_length_dev = 100


file_paternal = "paternal_alleles.txt"
file_maternal = "maternal_alleles.txt"
file_fetal = "fetal_alleles.txt"
file_samples = "plasma_samples.txt"

fpaternal = open(file_paternal, 'w')
fmaternal = open(file_maternal, 'w')
ffetal = open(file_fetal, 'w')
fsamples = open(file_samples, 'w')

#generate the haplotypes
M = []
P = []
for i in range(length):
    mnew = []
    pnew = []
    while mnew == pnew:
        alleles_set = []
        alleles_set.append(random.choice(nucleotides))
        alleles_set.append(random.choice(nucleotides))
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
  
inheritance_patterns = []
for mats in range(3):
    for pats in range(3):
        if mats+pats in [0,4]: continue
        inheritance_patterns.append([mats, pats]) 
        
#generate fetal alleles with CNVs
F = []
i = 0
X = 3
while i < length:
    #pick a type of the event (inheritence pattern)
    if not X == 3:
        X = random.choice([X, 3, 3])
    else:
        X = random.choice([0,1,2,3,4,5,6]) 
        
    #choose how many SNP the event will span
    #size = random.randrange(cnv_length_min, cnv_length_max)
    size = int(random.gauss(cnv_length_mu, cnv_length_dev))
        
    for j in range(size):
        if i == length: break
        fnew = []
        
        for k in range(inheritance_patterns[X][0]):
            fnew.append( random.choice(M[i]) )
        for k in range(inheritance_patterns[X][1]):
            fnew.append( random.choice(P[i]) )
        '''    
        if X==0: #ok
            fnew.append( random.choice(M[i]) )
            fnew.append( random.choice(P[i]) )
            
        elif X==1: #duplication
            fnew.append( random.choice(M[i]+P[i]) )
            fnew.append( random.choice(M[i]) )
            fnew.append( random.choice(P[i]) )
            
        elif X==2: #deletion
            fnew.append( random.choice(M[i]+P[i]) )
            
        elif X==3: #duplication+deletion
            if random.random()>0.5:
                fnew.append( random.choice(M[i]) )
                fnew.append( random.choice(M[i]) )
            else:
                fnew.append( random.choice(P[i]) )
                fnew.append( random.choice(P[i]) )
        '''        
        F.append(fnew)
        print >>ffetal, " ".join(fnew), X
        i += 1

#generate "reads"
for i in range(length):
    N = int(round(random.gauss(mean_coverage, dev_coverage)))
    samples = [0, 0, 0, 0]
    num_samples = N*mixture/2. * len(F[i]) + N*(1-mixture)/2. * len(M[i])
    for j in range(int(round(num_samples))):
        if random.random() <= 0.01:
            samples[random.choice(range(4))] += 1
        elif random.random() > (mixture * len(F[i])/2.) / (1.-mixture + mixture * len(F[i])/2.): #mixture/2. * len(F[i]): #from maternal alleles
            samples[nucleotides.index(random.choice(M[i]))] += 1
        else: #from fetal alleles
            samples[nucleotides.index(random.choice(F[i]))] += 1
    print >>fsamples, samples[0], samples[1], samples[2], samples[3]
    




