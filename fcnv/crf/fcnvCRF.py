import random
import math
#import numpypy as np
import itertools
import copy
import functools
from bisect import bisect_left

class memoized(object):
       """Decorator that caches a function's return value each time it is called.
       If called later with the same arguments, the cached value is returned, and
       not re-evaluated.
       """
       def __init__(self, func):
          self.func = func
          self.cache = {}
          
       def __call__(self, *args):
          try:
             return self.cache[args]
          except KeyError:
             value = self.func(*args)
             self.cache[args] = value
             return value
          except TypeError:
             # uncachable -- for instance, passing a list as an argument.
             # Better to not cache than to blow up entirely.
             print 'memo arg TypeError', args
             return self.func(*args)
             
       def __repr__(self):
          """Return the function's docstring."""
          return self.func.__doc__
          
       def __get__(self, obj, objtype):
          """Support instance methods."""
          fn = functools.partial(self.__call__, obj)
          fn.reset = self.reset
          return fn
          
       def reset(self):
          self.cache = {}

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
        '''
        dic = ["1", "2"]
        nms = ["M", "P"]
        result = ""
        for h in range(2):
            for x in self.phased_pattern[h]:
                result += " " + nms[h] + dic[x]
            if h == 0: result += " |"
        return result
        '''
        
    def __eq__(self, other):
        return self.inheritance_pattern == other.inheritance_pattern and \
            self.phased_pattern == other.phased_pattern

class FCNV(object):
    """Hidden Markov Model for fetal CNV calling"""
    
    nucleotides = ['A', 'C', 'G', 'T']
    
    def __init__(self, crfParams, positions, cnv_prior, use_prior):
        """Initialize new FCNV object"""
        super(FCNV, self).__init__()
        
        #store genomic positions of the SNPs
        self.positions = positions
        self.possible_distbins = range(5)
        if self.positions != None:
            self.distBin = [ self.getDistbin(x) for x in range(len(positions)+3)]
        
        #run intern tests
        self.neg_inf = float('-inf')
        self.inf = float('inf')
        self.nan = float('nan')
        self.testLogSum()
        
        #cache for log-likelihood values
        self.logLikelihoodCache = {}
        #self.distributionCache = {}
           
        '''#precompute lookup table for logSum function: 
        self.logTable1 = []
        for i in range(501):
            self.logTable1.append(math.log(1 + math.exp(-i/100.)) )
        self.logTable2 = []
        for i in range(51):
            self.logTable2.append(math.log(1 + math.exp(-i)) )
        '''
        
#        #hard-set hyper-parameters and default parameters
#         #weight vector for node features
#        self.unaryWeights = [1.]
#        self.unaryFeaturesList = [self.getUnaryF0]
#         #weight vector for edge features
#        #self.binaryWeights = [0.9799, 0.003, 0.0001, 0.9799, 0.001, 0.0001]
#        self.binaryWeights = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
#        self.binaryFeaturesList = [self.getBinaryF0, self.getBinaryF1, self.getBinaryF2, self.getBinaryF3, self.getBinaryF4, self.getBinaryF5]
#         #minimal epsilon weight for a feature
#        self.epsWeight = 0.0001
#         #training hyperparameters
#        self.sigmaSqr = 0.01
#        self.omega = 0.00005
        
        #generate inheritance patterns
        self.inheritance_patterns = []
        self.IPtoID = {}
        for mats in range(3):
            for pats in range(3):
                if mats+pats in [0,4]: continue
                self.IPtoID[(mats, pats)] = len(self.inheritance_patterns)
                self.inheritance_patterns.append((mats, pats))
        
               
        #generate phased variants of the inheritance patterns and list of HMM States
        self.states = []
        self.phased_patterns = []
        self.max_component_size = 0
        for ip in self.inheritance_patterns:
            num_phased_states = 0
            for mPhased in itertools.combinations_with_replacement([0,1], ip[0]):
                for pPhased in itertools.combinations_with_replacement([0,1], ip[1]):
                    phased_pattern = (tuple(mPhased), tuple(pPhased))
                    self.phased_patterns.append(phased_pattern)
                    new_state = HMMState(ip, phased_pattern)
                    self.states.append(new_state)
                    num_phased_states += 1
            self.max_component_size = max(self.max_component_size, num_phased_states)
        
#        for i, s in enumerate(self.states):
#            print i, s
        
        # CNV PRIOR
        #CNV per position prior
        self.use_prior = use_prior
        self.cnv_prior = cnv_prior
        print "Use priors: ", self.use_prior
        # biases
        self.cnv_prior_biases = [0, 0, 0]
        for state in self.states:
            copy_number = sum(state.inheritance_pattern) - 1
            self.cnv_prior_biases[copy_number] += 1
        
        #set hyper-parameters and CRF parameters
        self.setCRFparams(crfParams)
                
        #generate silent states
        self.states.append( HMMState("s", "s") ) #start state
        self.states.append( HMMState("t", "t") ) #end state

        #print the table
        num_states = self.getNumStates()
        #for k in xrange(num_states):
        #    for l in xrange(num_states):
        #        print "%.4f" % (math.exp(self.transitions[k][l])),
        #    print " "
        self.testEdgePotential()
        
        
    def binDist(self, dist):
        
        if 0 <= dist < 120:
            return 0
        elif dist < 320:
            return 1
        elif dist < 730:
            return 2
        elif dist < 3000:
            return 3
        else:
            return 4
    
    def getDistbin(self, pos):
        
        if 1 < pos < len(self.positions):
            return self.binDist(self.positions[pos-1] - self.positions[pos-2])
        else:
            return 0
            
    def setCRFparams(self, crfParams):
        """Set the given CRF parameters to global variables
        """
         #weight vector for node features
        self.unaryWeights = crfParams['unaryWeights']
        self.unaryFeaturesList = [self.getUnaryF0]
        if self.use_prior:
            self.unaryFeaturesList.append(self.getUnaryF1)
         #weight vector for edge features
        self.binaryWeights = crfParams['binaryWeights']
        self.binaryFeaturesList = []
        for distbin in self.possible_distbins:
            self.binaryFeaturesList.extend(self.makeBinary(distbin))
            
         #minimal epsilon weight for a feature
        self.epsWeight = crfParams['epsWeight']
         #training hyperparameters
        self.sigmaSqr = crfParams['sigmaSqr']
        self.omega = crfParams['omega']
        
    def encodeCRFparams(self):
        """Return dictionary with currently used CRF parameters
        """
        crfParams = {}
        crfParams['unaryWeights'] = self.unaryWeights
        crfParams['binaryWeights'] = self.binaryWeights
        crfParams['epsWeight'] = self.epsWeight
        crfParams['sigmaSqr'] = self.sigmaSqr
        crfParams['omega'] = self.omega
        return crfParams
            
    @memoized     
    def getNumIP(self):
        """Return number of inheritance patterns"""
        return len(self.inheritance_patterns)
    
    @memoized     
    def getNumPP(self):
        """Return number of phased inheritance patterns"""
        return len(self.phased_patterns)
    
    @memoized    
    def getNumStates(self):
        """Return number of HMM states per one SNP position"""
        return len(self.states)   
    
    @memoized 
    def isNormal(self, state):
        return state.inheritance_pattern == (1,1)
    
    @memoized   
    def isReal(self, state):
        return isinstance(state.inheritance_pattern, tuple) and isinstance(state.phased_pattern, tuple)
    
    def areRecombination(self, state1, state2):
        if state1.inheritance_pattern != state2.inheritance_pattern:
            return False
        eq = 0
        for hap in range(2):
            eq += (state1.phased_pattern[hap] == state2.phased_pattern[hap])
        return (eq > 0)
    
    @memoized
    def getCorrespondingInState(self, state):
        num_real_states = self.getNumPP()
        for st_id, st in enumerate(self.states[num_real_states:]):
            if st.inheritance_pattern == state.inheritance_pattern and \
               st.phased_pattern == "in":
                return num_real_states + st_id, st
        return None
    
    @memoized    
    def getCorrespondingOutState(self, state):
        num_real_states = self.getNumPP()
        for st_id, st in enumerate(self.states[num_real_states:]):
            if st.inheritance_pattern == state.inheritance_pattern and \
               st.phased_pattern == "out":
                return num_real_states + st_id, st
        return None
    
    @memoized    
    def getStartState(self):
        for i, state in enumerate(self.states):
            if state.inheritance_pattern == "s": return i, state
        return None
    
    @memoized   
    def getExitState(self):
        for i, state in enumerate(self.states):
            if state.inheritance_pattern == "t": return i, state
        return None
    
    def areAdjacent(self, state1, state2):
        """Find out whether the two states are adjacent.
        Adjacency is defined as hamming distance <= 1.
        
        """
        #convert from tuples to lists
        s = list(state1.phased_pattern)
        t = list(state2.phased_pattern)
        for i in range(2):
            s[i] = list(s[i])
            t[i] = list(t[i])
          
        #enumerate the  neighborhood of s with hamming distance 1
        #and look for discovery of t
        match = (s == t)
        for h in range(2): #over both haplotypes
            for i, x in enumerate(s[h]): #over all alleles in a haplotype
                temp = copy.deepcopy(s)
                temp[h][i] = (x + 1) % 2 #flip this allele
                match += (temp == t)
                
                del temp[h][i] #delete this allele
                match += (temp == t)
            for x in range(2):
                temp = copy.deepcopy(s)
                temp[h].append(x) #add a new allele
                match += (temp == t)
                
        return (match >= 1)
    
    def logNormalize(self, l):
        """Normalize the given list of log-probabilities.
        
        Normalizes the list in place and returns the used scale factor.
        
        """
        sum_ = self.neg_inf
        for x in l: sum_ = self.logSum(sum_, x)
        for i in range(len(l)): 
            if l[i] != self.neg_inf: l[i] -= sum_
        return sum_
        
    def logSum(self, x, y):
        """Return sum of two numbers in log space
        i.e. equivalent to log(exp(x)+exp(y))
        
        """
        '''
        def precise(x):
            return math.log(1 + math.exp(x) )

        def lookup(x):
            #return math.log(1 + math.exp(x) )
            x = -x
                    
            if x < 5:
                x *= 100
                fx = int(math.floor(x))
                return (x-fx)*(self.logTable1[fx+1] - self.logTable1[fx]) + self.logTable1[fx]
            elif x < 50:
                fx = int(math.floor(x))
                return (x-fx)*(self.logTable2[fx+1] - self.logTable2[fx]) + self.logTable2[fx]
            else: return 0.
        '''
        
        if x == self.neg_inf and y == self.neg_inf: return self.neg_inf
        elif x == self.inf and y == self.inf: return self.inf
        elif math.isnan(x) or math.isnan(y): return self.nan
        
        a = max(x, y)
        b = min(x, y)
        return a + math.log(1 + math.exp(b-a) )
    
    def testLogSum(self):
        # From the numpy documentation
        prob1 = math.log(1e-50)
        prob2 = math.log(2.5e-50)
        prob12 = self.logSum(prob1, prob2)
        assert prob12 == -113.87649168120691
        assert math.exp(prob12) == 3.5000000000000057e-50
        assert self.logSum(0, 0) == math.log(2)
        assert self.logSum(float('-inf'), 0) == 0
        #assert self.logSum(12345678, 12345678) == float('inf')

        assert math.isnan(self.logSum(float('nan'), 1))
        assert math.isnan(self.logSum(1, float('nan')))
        assert math.isnan(self.logSum(float('nan'), float('inf')))
        assert math.isnan(self.logSum(float('inf'), float('nan')))
        assert self.logSum(float('-inf'), float('-inf')) == float('-inf')
        assert self.logSum(float('-inf'), float('inf')) == float('inf')
        assert self.logSum(float('inf'), float('-inf')) == float('inf')
        assert self.logSum(float('inf'), float('inf')) == float('inf')
    
    #memoized
    def allelesDistribution(self, maternal, fetal, mix):
        """Compute nucleotides distribution in maternal plasma for a position with
        given maternal and fetal alleles, assuming the given fetal admixture.
        
        Arguments:
        maternal -- maternal alleles
        paternal -- paternal alleles
        mix -- fetal admixture ratio
        returns list(floats) -- list of nucleotides probabilities
        
        >>> f = FCNV()
        >>> f.allelesDistribution(['A', 'A'], ['T', 'A'], 0.1)
        [0.9230769230769231, 0.009615384615384616, 0.009615384615384616, 0.057692307692307696]
        """
        '''#Caching:
        code = (tuple(maternal), tuple(fetal), mix)
        val = self.distributionCache.get(code)
        if val != None: return val
        '''

        adjusted_fetal_admix = mix/2. * len(fetal)
        adjusted_maternal_admix = (1.-mix)/ 2. * len(maternal)
        cmix = adjusted_fetal_admix / (adjusted_maternal_admix + adjusted_fetal_admix)
        dist = []
        for nuc in self.nucleotides:
            p = maternal.count(nuc) / float(len(maternal)) * (1.-cmix)
            p += fetal.count(nuc) / float(len(fetal)) * (cmix)
            #if p < 0.01: p = 0.01
            dist.append(p)
            
        #normalize
        #summ = sum(dist)
        #dist = [dist[i] / summ for i in range(len(dist)) ]
        
        '''self.distributionCache[code] = dist'''
        return dist
    
    def allelesDistributionExtended(self, Malleles, Mcounts, Palleles, Pcounts, pattern, mix):
        """Compute nucleotides distribution in maternal plasma for a position with
        given maternal and fetal alleles, assuming the given fetal admixture.
        
        Arguments:
        maternal -- maternal alleles
        paternal -- paternal alleles
        mix -- fetal admixture ratio
        returns list(floats) -- list of nucleotides probabilities
        """

#        Falleles = []
#        #maternaly inherited
#        for mHpatt in pattern[0]:
#            Falleles.append(Malleles[mHpatt])
#        #paternaly inherited
#        for pHpatt in pattern[1]:
#            Falleles.append(Palleles[pHpatt])
            
        numFalleles = float(len(pattern[0]) + len(pattern[1]))
        numMinherited = len(pattern[0])
        numPinherited = len(pattern[1])
        numMalleles = float(len(Malleles))
        sumMcount = float(sum(Mcounts))
        sumPcount = float(sum(Pcounts))
        
        #adjust mixture to currently considered event (deletion, duplication, normal)
        adjusted_fetal_admix = mix/2. * numFalleles
        adjusted_maternal_admix = (1.-mix)/2. * numMalleles
        cmix = adjusted_fetal_admix / (adjusted_maternal_admix + adjusted_fetal_admix)
        
        dist = {}
        for nuc in self.nucleotides: dist[nuc] = 0.
        alphaM = 32/2. + sum(Mcounts)/2.
        alphaP = 39/2. + sum(Pcounts)/2.
        
        #fraction that is trully from mother DNA
        for i, nuc in enumerate(Malleles):
            dist[nuc] += (alphaM + Mcounts[i]) / float(2*alphaM + sumMcount) * (1.-cmix)
            
        #fraction that is fetal but maternaly inherited
        for mHpatt in pattern[0]:
            nuc = Malleles[mHpatt]
            if numMinherited == 2 and pattern[0][0] != pattern[0][1]:
                #if both haplotypes are inherited, use the maternal seq. ratio
                dist[nuc] += (alphaM + Mcounts[mHpatt]) / float(2*alphaM + sumMcount) * (cmix*2./numFalleles)
            else:
                dist[nuc] += cmix / numFalleles
            
        #fraction that is fetal but paternaly inherited
        for pHpatt in pattern[1]:
            nuc = Palleles[pHpatt]
            # dist[nuc] += cmix / numFalleles
            # use paternal ratio correction
            if numPinherited == 2 and pattern[1][0] != pattern[1][1]:
                #if both haplotypes are inherited, use the paternal seq. ratio
                dist[nuc] += (alphaP + Pcounts[pHpatt]) / float(2*alphaP + sumPcount) * (cmix*2./numFalleles)
            else:
                dist[nuc] += cmix / numFalleles
        
        dist_list = [dist[nuc] for nuc in self.nucleotides]
        #print pattern, Malleles, Mcounts, Palleles, Pcounts,':', dist_list, sum(dist_list)
        
        #normalize
        #summ = float(sum(dist_list))
        #dist_list = [dist_list[i] / summ for i in range(len(dist_list)) ]

        return dist_list
    
    def logGaussian(self, x, mus, cov_diagonal):
        '''
        log probability of x ~ N(mu, cov)
        '''
        D = len(x)
        N = sum(x)
        xm = [ float(x[i] - mus[i]) for i in range(D) ]
        
        det = 1.
        for i in range(D): det *= cov_diagonal[i]
            
        px = - D/2. * math.log(2.*math.pi)
        px += - math.log(det) /2.
        
        #compute the 'exponent' part #px += -(np.dot(np.dot(xm.T, inv), xm)) /2.
        part1 = [ xm[i] / cov_diagonal[i] for i in range(D) ]
        part2 = sum([ part1[i] * xm[i] for i in range(D) ])
        px += - part2 /2.
        
        return px
    
    def logDirMulti(self, rho, xs, ps):
        """
        log probability of x ~ BetaBinomial(xs, ps, rho)
        xs -- counts in individual categories
        ps -- means (probabilities) of individual categories; ps sum to 1
        rho -- overdispersion parameter of the BB distribution
        """
        
        def gammaln(n):
            """Compute logarithm of Euler's gamma function for discrete values."""
            if n < 1:
                return float('inf')
            if n < 3:
                return 0.0
            c = [76.18009172947146, -86.50532032941677, \
                 24.01409824083091, -1.231739572450155, \
                 0.001208650973866179, -0.5395239384953 * 0.00001]
            x, y = float(n), float(n)
            tm = x + 5.5
            tm -= (x + 0.5) * math.log(tm)
            se = 1.0000000000000190015
            for j in range(6):
                y += 1.0
                se += c[j] / y
            return -tm + math.log(2.5066282746310005 * se / x)
        
        def logFactorial(x):
            """Calculate ln(x!).
                
            Arguments:
            x -- list(floats)
            returns list(floats)
                
            """
            if isinstance(x, tuple):
                res = []
                for val in x:
                    res.append(gammaln(val+1))
                return tuple(res)
            else: 
                return gammaln(x+1)
        
        res = 0.
        for i in range(len(xs)):
            for r in range(1, xs[i]+1):
                res += math.log( ps[i]*(1.-rho) + (r-1.)*rho)
        for r in range(1, sum(xs)+1):
            res -= math.log( (1.-rho) + (r-1.)*rho )
            
        res += logFactorial(sum(xs)) - sum(logFactorial(xs))
        return res
    
    def getAlleleCounts(self, pos_ind, nuc_counts, maternal_alleles, paternal_alleles,\
           maternal_sq_counts, paternal_sq_counts, mix, state_id, parameterStats):
        
        num_real_states = self.getNumPP()      
        state = self.states[state_id]

        ps = self.allelesDistributionExtended(maternal_alleles, maternal_sq_counts, paternal_alleles, paternal_sq_counts, state.phased_pattern, mix)
        
        #identify the two possible alleles
        expinherited = [0] * 4
        for x in maternal_alleles + paternal_alleles:
            expinherited[self.nucleotides.index(x)] = 1
        
        seen = [0] * 4
        for x in range(4):
            if nuc_counts[x] != 0: seen[x] = 1
            
        for x in range(4):
            seen[x] = int(seen[x] or expinherited[x])
        
        als = [i for i in range(4) if seen[i] == 1]
        
        #print nuc_counts, '|', als, '|', seen
        if len(als) > 2: 
            print "WARNING: trialelic pos"
            return
        if len(als) == 1: #homozygous position, so just use whatever second allele (has 0 count)
            for x in range(4):
                if x != als[0]: 
                    als.append(x)
                    break

        
        N = float(sum(nuc_counts))
        m = len(als)
        dmps = [0.] * m
        for i in range(m):
            dmps[i] = ps[als[i]]
        dmps = [ x/sum(dmps) for x in dmps]
        
        #make the one with higher p to be the first one
        if dmps[1] > dmps[0]:
            dmps[0], dmps[1] = dmps[1], dmps[0]
            als[0], als[1] = als[1], als[0]
            
        expP = dmps[0]
        alleleAcount = nuc_counts[als[0]]
        alleleBcount = nuc_counts[als[1]]
            
        if expP not in parameterStats:
            parameterStats[expP] = list()
        parameterStats[expP].append( (alleleAcount, alleleBcount) )
        
        #print expP, (alleleAcount, alleleBcount), "|", nuc_counts, ps
    
    def logLHGivenStateWCoverage(self, pos_ind, nuc_counts, maternal_alleles, paternal_alleles,\
           maternal_sq_counts, paternal_sq_counts, mix, state):
        
        if state.inheritance_pattern == (2, 0) or state.inheritance_pattern == (0, 2): return -50

        
        pattern = state.phased_pattern
        ps = self.allelesDistributionExtended(maternal_alleles, maternal_sq_counts, paternal_alleles, paternal_sq_counts, pattern, mix)

#        #use Mixture of Independent Gaussians
#        N = sum(nuc_counts)
#        avg_doc = 70. #TODO: make this a parameter
#        norm_coef = avg_doc / N
#        
#        mus = [ avg_doc * ps[i] for i in range(4) ]
#        cov_diagonal = [ max(0.8, (N * ps[x]) * norm_coef**2) for x in range(4)]
#        nuc_counts = [ x * norm_coef for x in nuc_counts ]
#        
#        ratios_prob = self.logGaussian(nuc_counts, mus, cov_diagonal)
#        result = ratios_prob

        #use Beta-Binomial distrib
        expinherited = [0] * 4
        for x in maternal_alleles + paternal_alleles:
            expinherited[self.nucleotides.index(x)] = 1
        #als = [i for i in range(4) if expinherited[i] == 1]
        
        seen = [0] * 4
        for x in range(4):
            if nuc_counts[x] != 0: seen[x] = 1
            
        for x in range(4):
            seen[x] = int(seen[x] or expinherited[x])
        
        als = [i for i in range(4) if seen[i] == 1]
        
        #print nuc_counts, '|', als, '|', seen
        if len(als) > 2: return -3
        if len(als) == 1:
            for x in range(4):
                if x != als[0]: 
                    als.append(x)
                    break
        
        m = len(als)
        eps = 0.01
        dmps = [0.] * m
        addEps = False
        for i in range(m):
            if ps[als[i]] < 0.0001: addEps = True
        for i in range(m):
            dmps[i] = ps[als[i]] + eps*int(addEps)
        dmps = [ x/sum(dmps) for x in dmps]
        
        #make the one with higher p to be the first one
        if dmps[1] > dmps[0]:
            dmps[0], dmps[1] = dmps[1], dmps[0]
            als[0], als[1] = als[1], als[0]
        
        rho = 0.0267  #0.01608  #estimated from train set by MoM estimator
        nc = [nuc_counts[als[i]] for i in range(m)]
        result = self.logDirMulti(rho, tuple(nc), tuple(dmps))
        
        #if result < -15: result = -15
        
        return result
    
    def estimateMixture(self, samples, M, P):
        """Estimate mixture of fetal genome in the plasma samples.
        
        Arguments:
        samples -- number of nucleotides sampled per individual SNP positions
        M -- maternal alleles
        P -- parental alleles
        returns float
        
        """
        mix = 0.
        num = 0
        estimates = []
        for i in xrange(len(samples)):
            #if mother is homozygous A and father has at least one alternative B
            if M[i][0] == M[i][1] and (P[i][0] != M[i][0] and P[i][1] != M[i][0]):
                #print M[i], P[i], samples[i]            
                type1 = M[i][0]
                type2 = P[i][0]
                if type1 == type2: type2 = P[i][1]
                try:
                    num_type2 = samples[i][self.nucleotides.index(type2)]
                except: 
                    print "ERROR unexpected allele: ", type2
                    
                #only if it looks like the child inherited the allele B
                if num_type2 != 0: 
                    sum_all = sum(samples[i])
                    num += 1
                    mix += (2.*num_type2)/sum_all
                    estimates.append((2.*num_type2)/sum_all)
                else:
                    num += 1.
                    estimates.append(0.)

        estimates.sort()
        return mix/num, estimates[len(estimates)/2], len(estimates)
        
    def estimateCoverage(self, samples):
        """Estimate coverage of the SNP loci in the plasma samples as empirical mean.
        """
        result = sum(map(sum, samples)) / float(len(samples))
        return result
    
    def score(self, labels, samples, M, P, MSC, PSC, mixture):
        score = 0
        for pos, sample in enumerate(samples):

            state = labels[pos]
            
            for w, f in self.unaryWeightFeaturesList:
                score += f(pos+1, samples, M, P, MSC, PSC, mixture, y)*w

            # transitions            
            if pos > 0:
                prev_state = labels[pos-1]
                score += self.binaryFeature(prev_state, state, pos)
            
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    def computeLLandGradient(self, labels, samples, M, P, MSC, PSC, mixture):
        """
        Compute the parameters likelihood, corresponding gradient, and update the weights
        """
        distBin = self.distBin
        #get forward and backward DP tables and log of p(X) -- the Z function
        fwd, pX1, scale = self.computeForward(samples, M, P, MSC, PSC, mixture)
        bck, pX2 = self.computeBackward(samples, M, P, MSC, PSC, mixture, scale)
        
        logZ = pX1 + sum(scale)
        
        if abs(pX1 - pX2) > 10e-8:
            print "p(X) from forward and backward DP doesn't match", pX1, pX2, sum(scale)
            return
        
        num_real_states = self.getNumPP()
        n = len(samples)
        
        #compute node marginals
        nodeMarginals = [ [0.]*num_real_states for y in range(n+1)]
        #logNodeMarginals = [ [self.neg_inf]*num_real_states for y in range(n+1)]
        for pos in xrange(1, n+1):
            #logSumm = self.neg_inf
            for s_id in xrange(num_real_states):
                logpStateAtPos = fwd[pos][s_id] + bck[pos][s_id] - pX1
                nodeMarginals[pos][s_id] = math.exp(logpStateAtPos)
                #logNodeMarginals[pos][s_id] = logpStateAtPos
                #logSumm = self.logSum(fwd[pos][s_id] + bck[pos][s_id] - pX1, logSumm)
            #print logSumm, math.exp(logSumm)
        
        #compute egde marginals
        edgePotential = self.getEdgePotential(self.binaryWeights)
        edgeMarginals = [ [[0.]*num_real_states for x in range(num_real_states)] for y in range(n+1)]
        
        for pos in xrange(1, n):
            nodePotential = self.getNodePotential(pos+1, self.unaryWeights, samples, M, P, MSC, PSC, mixture)
            logSumm = self.neg_inf
            for s1_id in xrange(num_real_states):
                for s2_id in xrange(num_real_states):
                    pEdgeAtPos = fwd[pos][s1_id] + edgePotential[s1_id][s2_id][distBin[pos+1]] + nodePotential[s2_id] + bck[pos+1][s2_id]
                    edgeMarginals[pos][s1_id][s2_id] = pEdgeAtPos
                    logSumm = self.logSum(pEdgeAtPos, logSumm)
            for s1_id in xrange(num_real_states):
                for s2_id in xrange(num_real_states):
                    edgeMarginals[pos][s1_id][s2_id] = math.exp(edgeMarginals[pos][s1_id][s2_id] - logSumm)
        
        #compute the likelihood of current parameters (weight vectors)           
        logLikelihood, = self.getScore(samples, M, P, MSC, PSC, mixture, labels)
        #for i, x in enumerate(labels):
        #    print labels[i], x, self.IPtoID[self.states[x].inheritance_pattern]
            
        print logLikelihood, logZ
        logLikelihood -= logZ
        
        wSqr = sum([ x*x  for x in self.unaryWeights + self.binaryWeights ])
        logLikelihood -= wSqr/(2.*self.sigmaSqr) #regularizator
        
        print "loglikelihood before param adjustment:", logLikelihood
        
        #compute gradients w.r.t. indiviudal features and update the current weights
         #unary features
        for i, f in enumerate(self.unaryFeaturesList):
            grad = 0.
            for pos in range(len(labels)):
                grad += math.exp(f(pos+1, samples, M, P, MSC, PSC, mixture, self.states[labels[pos]]))
                expect = 0.  #self.neg_inf
                for s_id in xrange(num_real_states):
                    expect += nodeMarginals[pos+1][s_id] * math.exp(f(pos+1, samples, M, P, MSC, PSC, mixture, self.states[s_id]))
                grad -= expect
            grad -= self.unaryWeights[i]/self.sigmaSqr #regularizator
            #update the current weights
            self.unaryWeights[i] += self.omega * grad
            self.unaryWeights[i] = max(self.unaryWeights[i], self.epsWeight)
            print "unary", i, self.unaryWeights[i]
            
         #binary features
        for i, f in enumerate(self.binaryFeaturesList):
            grad = 0.
            for pos in xrange(len(labels)-1):
                grad += f(labels[pos], self.states[labels[pos]], labels[pos+1], self.states[labels[pos+1]], distBin[pos+2])
                expect = 0.
                for s1_id in xrange(num_real_states):
                    for s2_id in xrange(num_real_states):
                        expect += edgeMarginals[pos+1][s1_id][s2_id] * f(s1_id, self.states[s1_id], s2_id, self.states[s2_id], distBin[pos+2])
                grad -= expect
            grad -= self.binaryWeights[i]/self.sigmaSqr #regularizator
            #update the current weights
            self.binaryWeights[i] += self.omega * grad
            self.binaryWeights[i] = max(self.binaryWeights[i], self.epsWeight)
            print "binary", i, self.binaryWeights[i]    
        
        return logLikelihood, self.encodeCRFparams()
       
    def generateWeightMatrixForMCC(self):
        """
        Creates 'num_real_states' x 'num_real_states' matrix of weights for the 
        individual correlation coefficients used in the MCC computation
        """
        num_real_states = self.getNumPP()  
        w = [ [0.] * num_real_states for x in range(num_real_states)]
        
        for s1_id, s1 in enumerate(self.states[:num_real_states]):
            for s2_id, s2 in enumerate(self.states[:num_real_states]):
                #if ground truth is a normal state
                if s1.inheritance_pattern == (1,1):
                    #the same state
                    if s1_id == s2_id:
                        w[s1_id][s2_id] = 0.
                    #recombination
                    elif s1.inheritance_pattern == s2.inheritance_pattern:
                        w[s1_id][s2_id] = 0.
                    #other inheritance pattern
                    else:
                        w[s1_id][s2_id] = 1.
                #else if ground truth is a CNV state
                else:
                    #the same state
                    if s1_id == s2_id:
                        w[s1_id][s2_id] = 1.
                    #recombination
                    elif s1.inheritance_pattern == s2.inheritance_pattern:
                        w[s1_id][s2_id] = 0.5
                    #other inheritance pattern
                    else:
                        w[s1_id][s2_id] = 1.
        
#        for i in range(len(w)):
#            for j in range(len(w[0])):
#                print w[i][j],
#            print ''
        
        return w

    def generateWeightMatrixForSWE(self):
        """
        Creates 'num_real_states' x 'num_real_states' matrix of weights for the 
        individual correlation coefficients used in the MCC computation
        """
        num_real_states = self.getNumPP()  
        w = [ [0.] * num_real_states for x in range(num_real_states)]
        
        for s1_id, s1 in enumerate(self.states[:num_real_states]):
            for s2_id, s2 in enumerate(self.states[:num_real_states]):

                #the same state
                if s1_id == s2_id:
                    w[s1_id][s2_id] = 0.0

                else:
                    #if ground truth is a normal state
                    if s1.inheritance_pattern == (1,1):

                        #recombination
                        if s1.inheritance_pattern == s2.inheritance_pattern:
                            w[s1_id][s2_id] = 0.0

                        #other inheritance pattern                        
                        else:
                            w[s1_id][s2_id] = 0.5

                    #else if ground truth is a CNV state
                    else:
                        
                        #recombination
                        if s1.inheritance_pattern == s2.inheritance_pattern:
                            w[s1_id][s2_id] = 0.5

                        #other inheritance pattern
                        else:
                            w[s1_id][s2_id] = 1.
                
        return w
        
        
    def getConfusionMatrix(self, labels, predicted_labels):
        """
        Compute multiclass confusion matrix
        """
        num_states = self.getNumPP()        
        confusionMatrix = [[0. for x in range(num_states)] for y in range(num_states)]
        
        for i in xrange(len(labels)):
            confusionMatrix[labels[i]][predicted_labels[i]] += 1
        
        return confusionMatrix
        
    def matthewsCorrelationCoefficientLoss(self, labels, predicted_labels):
        """
        Compute the Matthews Correlation Coefficient from the Confusion Matrix of
        labels and predicted_labels. Returns a float.
        """
        
#        print "MCC:"
#        print len(labels), len(predicted_labels)
#        ok_pred = 0
#        for i in range(len(labels)):
#            print i, '\t', labels[i], predicted_labels[i]
#            if labels[i] == predicted_labels[i]: ok_pred +=1
#        print "MCC END"
#        print "ok pred: ", ok_pred
        
        num_states = self.getNumPP()
        confusionMatrix = self.getConfusionMatrix(labels, predicted_labels)
        covWeightMatrix = self.generateWeightMatrixForMCC()
                
        print "Confusion Matrix:"
        for x in range(len(confusionMatrix)):
            for y in range(len(confusionMatrix[x])):
                print int(confusionMatrix[x][y]),
            print ''
        
        #reweight the confusion matrix
        for i in range(len(confusionMatrix)):
            for j in range(len(confusionMatrix[i])):
                confusionMatrix[i][j] *= covWeightMatrix[i][j]
        
        print "Confusion Matrix AFTER REWEIGHTING:"
        for x in range(len(confusionMatrix)):
            for y in range(len(confusionMatrix[x])):
                print int(confusionMatrix[x][y]),
            print ''
        
        cov = 0.
        varx = 0.
        vary = 0.
        
        for k in xrange(num_states):
            varxfact1 = 0.0
            varxfact2 = 0.0
            varyfact1 = 0.0
            varyfact2 = 0.0
            
            for l in xrange(num_states):                
                varxfact1 += confusionMatrix[l][k]
                varyfact1 += confusionMatrix[k][l]
                
                for m in xrange(num_states):
                    cov += confusionMatrix[k][k] * confusionMatrix[m][l] - confusionMatrix[l][k]*confusionMatrix[k][m]
                    
                    if l != k:
                        varxfact2 += confusionMatrix[m][l]
                        varyfact2 += confusionMatrix[l][m]
                    
            varx += varxfact1*varxfact2
            vary += varyfact1*varyfact2
            
        try:
            mcc = cov/(math.sqrt(varx)*math.sqrt(vary))
            print "mcc: ", mcc
        except Exception as e:
            print "Error in MCC computation:", e, cov, varx, vary
            mcc = 0.0
            
        return 1- mcc

    def sumWeightedErrors(self, labels, predicted_labels):
        """
        Computes the dot product between the confusion matrix
        and the matrix returned by generateWeightMatrixForSWE.
        """
        
        
        num_states = self.getNumPP()
        confusionMatrix = self.getConfusionMatrix(labels, predicted_labels)
        w = self.generateWeightMatrixForSWE()
        N, M = len(w), len(w[0])
        return sum(sum(w[i][j]*confusionMatrix[i][j] for j in range(M)) for i in range(N))

    def confusionEntropy(self, labels, predicted_labels):
        """
        Compute the Matthews Correlation Coefficient from the Confusion Matrix of
        labels and predicted_labels. Returns a float.
        """
        
        
        confusionMatrix = self.getConfusionMatrix(labels, predicted_labels)
        
        CEN = 0.
        num_states = self.getNumPP()
        
        for j in xrange(num_states):
            pj = 0.0
            pj_normalization = 0.0

            pjk_j = 0.0
            pjk_j = 0.0
            
            for k in xrange(num_states):
                pj += confusionMatrix[j][k] + confusionMatrix[k][j]
                pj_normalization += confusionMatrix[k][l]
                
            e = 0.
            
            for k in xrange(num_states):
                pkj_j = confusionMatrix[k][j]/pj
                pjk_j = confusionMatrix[j][k]/pj
                
                if k != j:
                    if pjk_j > 0.0:
                        e += - pjk_j * math.log(pjk_j, 2*num_states - 2) 
                    if pkj_j > 0.0:
                        e += - pkj_j * math.log(pkj_j, 2*num_states - 2)
                
            pj = 0.5*pj/pj_normalization
            CEN += pj*e
            
        return CEN
    
    def getScore(self, samples, M, P, MSC, PSC, mixture, *labelslists):
        """
        Computes logScore of (multiple) lists of sequence labeles
        """
        distBin = self.distBin
        edgePotential = self.getEdgePotential(self.binaryWeights)
        
        #scores = [self.neg_inf] * len(labelslists)
        scores = [0.] * len(labelslists)
        length = len(labelslists[0])
        for pos in xrange(length):
            # print labels[pos], predicted_labels[pos]
            #real positions are <1..n+1)
            nodePotential = self.getNodePotential(pos+1, self.unaryWeights, samples, M, P, MSC, PSC, mixture)
            
            for i in range(len(scores)):
                #pairwise factor value
                if pos+1 < length:
                    scores[i] += edgePotential[labelslists[i][pos]][labelslists[i][pos+1]][ distBin[pos+2] ]
                #unary factor value
                scores[i] += nodePotential[labelslists[i][pos]]
                #if nodePotential[labelslists[i][pos]] == self.neg_inf: print pos, ":", labelslists[i][pos]
                
        #for i in range(len(scores)):
        #    scores[i] = math.exp(scores[i])
            
        return scores

    def computeCountsTable(self, labels, samples, M, P, MSC, PSC, mixture, parameterStats):
        n = len(samples)
        
        # compute stats of the main distrib. parameter        
        for pos in xrange(1, n+1):
            #real positions are <1..n+1), pos 1 is the base case    
            self.getAlleleCounts(pos-1, samples[pos-1], M[pos-1], P[pos-1], MSC[pos-1], PSC[pos-1], mixture, labels[pos-1], parameterStats)
        
        return parameterStats
        
    def computeLLandMaxMarginUpdate(self, labels, samples, M, P, MSC, PSC, mixture, C, compute_postloss=True):
        """
        Compute the parameters likelihood, max margin update, and update the weights
        """
        n = len(samples)
        
                
        # SCORE OF THE TRUE LABELS
        predicted_labels, _ = self.viterbiPath(samples, M, P, MSC, PSC, mixture, inferHighLevelLabels=False)
        groundtruth_score, prediction_score = self.getScore(samples, M, P, MSC, PSC, mixture, labels, predicted_labels)
        
        print predicted_labels[:100]
        print labels[:100]
        print "labels scores:", groundtruth_score, prediction_score
        
        # groundtruth_score = self.neg_inf
        # prediction_score = self.neg_inf
        # 
        # for pos in xrange(len(labels)):
        #     # print labels[pos], predicted_labels[pos]
        #     #pairwise factor value
        #     if pos+1 < len(labels):
        #         groundtruth_score = self.logSum(groundtruth_score, edgePotential[labels[pos]][labels[pos+1]])
        #         prediction_score = self.logSum(prediction_score, edgePotential[predicted_labels[pos]][predicted_labels[pos+1]])
        #         
        #     #real positions are <1..n+1)
        #     nodePotential = self.getNodePotential(pos+1, self.unaryWeights, samples, M, P, MSC, PSC, mixture)
        #     groundtruth_score = self.logSum(groundtruth_score, nodePotential[labels[pos]])
        #     prediction_score = self.logSum(prediction_score, nodePotential[predicted_labels[pos]])
        #             
        # groundtruth_score = math.exp(groundtruth_score)
        # prediction_score = math.exp(prediction_score)
            
        # COMPUTE FEATURES
        groundtruthUnaryFeatures = self.getUnaryFeatures(labels, samples, M, P, MSC, PSC, mixture) 
        groundtruthBinaryFeatures = self.getBinaryFeatures(labels, samples, M, P, MSC, PSC, mixture) 
        predictionUnaryFeatures = self.getUnaryFeatures(predicted_labels, samples, M, P, MSC, PSC, mixture) 
        predictionBinaryFeatures = self.getBinaryFeatures(predicted_labels, samples, M, P, MSC, PSC, mixture)

        # COMPUTE SQUARED DISTANCE
        squared_feature_distance = 0
        for i in xrange(len(predictionUnaryFeatures)):
            squared_feature_distance += (groundtruthUnaryFeatures[i] - predictionUnaryFeatures[i]) ** 2

        print "sqr unary feature dist: ", squared_feature_distance
        
        for i in xrange(len(predictionBinaryFeatures)):
            squared_feature_distance += (groundtruthBinaryFeatures[i] - predictionBinaryFeatures[i]) ** 2
        
        print "total sqrt f.dist.:", squared_feature_distance
        
        # COMPUTE LOSS
        loss = self.matthewsCorrelationCoefficientLoss(labels, predicted_labels)
        
        # COMPUTE TAU
        print "loss: {0}".format(loss)
        print "prediction_score: {0}".format(prediction_score)
        print "groundtruth_score: {0}".format(groundtruth_score)
        print "squared_feature_distance: {0}".format(squared_feature_distance)
        
        tau = min(C, (prediction_score - groundtruth_score + loss)/squared_feature_distance)
        
        print "tau: {0}".format(tau)
        
        for i, f in enumerate(self.unaryFeaturesList):
            update = tau * (groundtruthUnaryFeatures[i] - predictionUnaryFeatures[i])
            self.unaryWeights[i] += update
            self.unaryWeights[i] = max(self.unaryWeights[i], self.epsWeight)
            
        for i, f in enumerate(self.binaryFeaturesList):
            update = tau * (groundtruthBinaryFeatures[i] - predictionBinaryFeatures[i])
            self.binaryWeights[i] += update
            self.binaryWeights[i] = max(self.binaryWeights[i], self.epsWeight)
            
        if compute_postloss:
            predicted_labels, _ = self.viterbiPath(samples, M, P, MSC, PSC, mixture, inferHighLevelLabels=False)
            postloss = self.matthewsCorrelationCoefficientLoss(labels, predicted_labels)
            postgroundtruth_score, postprediction_score = self.getScore(samples, M, P, MSC, PSC, mixture, labels, predicted_labels)
        else:
            postloss, postgroundtruth_score, postprediction_score = None, None, None
            
        return groundtruth_score, prediction_score, loss, self.encodeCRFparams(), postgroundtruth_score, postprediction_score, postloss
    
    def getUnaryFeatures(self, labels, samples, M, P, MSC, PSC, mixture):
        
        unaryFeatures = []
        for i, f in enumerate(self.unaryFeaturesList):
            feature = self.neg_inf
            for pos in xrange(len(labels)):
                feature = self.logSum(feature, f(pos+1, samples, M, P, MSC, PSC, mixture, self.states[labels[pos]]))
            unaryFeatures.append(math.exp(feature))
            
        return unaryFeatures

    def getBinaryFeatures(self, labels, samples, M, P, MSC, PSC, mixture):

        binaryFeatures = []
        for i, f in enumerate(self.binaryFeaturesList):
            feature = 0.
            for pos in xrange(len(labels)-1):
                feature += f(labels[pos], self.states[labels[pos]], labels[pos+1], self.states[labels[pos+1]], self.distBin[pos+2])
            binaryFeatures.append(feature)
        return binaryFeatures
    
    def getUnaryF0(self, pos, samples, M, P, MSC, PSC, mixture, state):
        #emission probability in the given state
        emis_logp = self.logLHGivenStateWCoverage(pos-1, samples[pos-1], M[pos-1], P[pos-1], MSC[pos-1], PSC[pos-1], mixture, state)
        return emis_logp

    def getUnaryF1(self, pos, samples, M, P, MSC, PSC, mixture, state):
        #copy number prior for the given state
        copy_number = sum(state.inheritance_pattern) - 1
        res = self.cnv_prior[pos-1][copy_number] - math.log(self.cnv_prior_biases[copy_number])
        return res
        
    def getNodePotential(self, pos, unaryWeights, samples, M, P, MSC, PSC, mixture):
        """
        Compute overall node potential for position @pos in the observation sequence
        
        return an array of 'num states' length representing the potential for 
        each possible label (state), by pooling all node features together
        given weights vector @unaryWeights
        """
        
        featuresList = self.unaryFeaturesList
        num_real_states = self.getNumPP()
        nodePot = [self.neg_inf] * num_real_states
        for state_id, state in enumerate(self.states[:num_real_states]):
            #sum log values of all unary features for this position and state/lable
            for i, f in enumerate(featuresList):
                featureValue = f(pos, samples, M, P, MSC, PSC, mixture, state)
                featureValue += math.log(unaryWeights[i])
                nodePot[state_id] = self.logSum(featureValue, nodePot[state_id])
            nodePot[state_id] = math.exp(nodePot[state_id])
            
        return nodePot

    def makeBinary(self, distbin):
        
        def getBinaryF0(sid1, state1, sid2, state2, curr_distbin):
            #for Normal states: pairwise energy of staying in the state
            if distbin == curr_distbin:
                if state1.inheritance_pattern == (1, 1):
                    if sid1 == sid2:
                        return 1
            return 0
        
        def getBinaryF1(sid1, state1, sid2, state2, curr_distbin):
            #for Normal states: pairwise energy of recombination of the same IP
            if distbin == curr_distbin:
                if state1.inheritance_pattern == (1, 1):
                    if sid1 != sid2 and state1.inheritance_pattern == state2.inheritance_pattern:
                        return 1     
            return 0
    
        def getBinaryF2(sid1, state1, sid2, state2, curr_distbin):
            #for Normal states: pairwise energy of going to other real states -- to CNV nodes
            if distbin == curr_distbin:
                if state1.inheritance_pattern == (1, 1):
                    if state1.inheritance_pattern != state2.inheritance_pattern and self.isReal(state2):
                        if state2.inheritance_pattern != (0, 2) and state2.inheritance_pattern != (2, 0):
                            return 1
            return 0
        
        def getBinaryF3(sid1, state1, sid2, state2, curr_distbin):
            #for CNV states: pairwise energy of staying in the state
            if distbin == curr_distbin:
                if state1.inheritance_pattern != (1, 1) and self.isReal(state1):
                    if sid1 == sid2:
                        return 1
            return 0
        
        def getBinaryF4(sid1, state1, sid2, state2, curr_distbin):
            #for CNV states: pairwise energy of recombination of the same IP
            if distbin == curr_distbin:
                if state1.inheritance_pattern != (1, 1) and self.isReal(state1):
                    if sid1 != sid2 and state1.inheritance_pattern == state2.inheritance_pattern:
                        result = 1
                        if state1.inheritance_pattern == (2, 1):
                            if state1.phased_pattern[0][0] ^ state1.phased_pattern[0][1] != state2.phased_pattern[0][0] ^ state2.phased_pattern[0][1]:
                                result = 0
                        elif state1.inheritance_pattern == (1, 2):
                            if state1.phased_pattern[1][0] ^ state1.phased_pattern[1][1] != state2.phased_pattern[1][0] ^ state2.phased_pattern[1][1]:
                                result = 0
                        return result
            return 0
    
        def getBinaryF5(sid1, state1, sid2, state2, curr_distbin):
            #for CNV states: pairwise energy of going to Normal states
            if distbin == curr_distbin:
                if state1.inheritance_pattern != (1, 1) and self.isReal(state1):
                    if state2.inheritance_pattern == (1, 1):
                        return 1
            return 0
        
        def getBinaryF6(sid1, state1, sid2, state2, curr_distbin):
            #for CNV states: pairwise energy of going to too diffrent recombination
            if distbin != curr_distbin: return 0
            if state1.inheritance_pattern != (1, 1) and self.isReal(state1):
                if sid1 != sid2 and state1.inheritance_pattern == state2.inheritance_pattern:
                    if state1.inheritance_pattern == (2, 1):
                        if state1.phased_pattern[0][0] ^ state1.phased_pattern[0][1] != state2.phased_pattern[0][0] ^ state2.phased_pattern[0][1]:
                            return 1
                    elif state1.inheritance_pattern == (1, 2):
                        if state1.phased_pattern[1][0] ^ state1.phased_pattern[1][1] != state2.phased_pattern[1][0] ^ state2.phased_pattern[1][1]:
                            return 1
            return 0
        
        return [getBinaryF0, getBinaryF1, getBinaryF2, getBinaryF3, getBinaryF4, getBinaryF5, getBinaryF6]
        
    def getEdgePotential(self, binaryWeights):
        """
        Compute overall edge potential shared across all edges
        
        return 'num states' x 'num states' x 'num DistBins' matrix representing the potential that
        pools all edge features together given weights vector @binaryWeights
        """
        num_real_states = self.getNumPP()
        num_states = self.getNumStates()
        
        featuresList = self.binaryFeaturesList
        edgePot = [[[self.neg_inf for db in self.possible_distbins] for x in range(num_states)] for y in range(num_states)]
        
        for dbid, db in enumerate(self.possible_distbins):
            for sid1, state1 in enumerate(self.states[:num_real_states]):
                for sid2, state2 in enumerate(self.states[:num_real_states]):
                    
                    # sum log values of all binary features 
                    # for pair of states/lables
                    for i, f in enumerate(featuresList):
                        featureValue = f(sid1, state1, sid2, state2, db)
                        if featureValue != 0 and edgePot[sid1][sid2][dbid] == self.neg_inf: edgePot[sid1][sid2][dbid] = 0.
                        
                        featureValue *= binaryWeights[i]
                        edgePot[sid1][sid2][dbid] += featureValue
            
                #constant to the exit state
                sid2 = self.getExitState()[0]
                edgePot[sid1][sid2][dbid] = 1.  
        
            #now generate constant pairwise energy from the start node
            for sid1, state1 in enumerate(self.states[num_real_states:]):
                sid1 += num_real_states
                #if it is the start node
                if state1.phased_pattern == "s":
                    for sid2, state2 in enumerate(self.states[:num_real_states]):
                        if state2.inheritance_pattern != (0, 2) and state2.inheritance_pattern != (2, 0):
                            edgePot[sid1][sid2][dbid] = 1.
                
        return edgePot

    def getEdgePotential(self, binaryWeights):
        """
        Compute overall edge potential shared across all edges
        
        return 'num states' x 'num states' x 'num DistBins' matrix representing the potential that
        pools all edge features together given weights vector @binaryWeights
        """
        num_real_states = self.getNumPP()
        num_states = self.getNumStates()
        
        featuresList = self.binaryFeaturesList
        edgePot = [[[self.neg_inf for db in self.possible_distbins] for x in range(num_states)] for y in range(num_states)]
        
        for dbid, db in enumerate(self.possible_distbins):
            for sid1, state1 in enumerate(self.states[:num_real_states]):
                for sid2, state2 in enumerate(self.states[:num_real_states]):
                    
                    # sum log values of all binary features 
                    # for pair of states/lables
                    for i, f in enumerate(featuresList):
                        featureValue = f(sid1, state1, sid2, state2, db)
                        if featureValue != 0 and edgePot[sid1][sid2][dbid] == self.neg_inf: edgePot[sid1][sid2][dbid] = 0.
                        
                        featureValue *= binaryWeights[i]
                        edgePot[sid1][sid2][dbid] += featureValue
            
                #constant to the exit state
                sid2 = self.getExitState()[0]
                edgePot[sid1][sid2][dbid] = 1.  
        
            #now generate constant pairwise energy from the start node
            for sid1, state1 in enumerate(self.states[num_real_states:]):
                sid1 += num_real_states
                #if it is the start node
                if state1.phased_pattern == "s":
                    for sid2, state2 in enumerate(self.states[:num_real_states]):
                        if state2.inheritance_pattern != (0, 2) and state2.inheritance_pattern != (2, 0):
                            edgePot[sid1][sid2][dbid] = 1.
                
        return edgePot
            
    def testEdgePotential(self):
        """
        Compute overall edge potential shared across all edges
        
        return 'num states' x 'num states' x 'num DistBins' matrix representing the potential that
        pools all edge features together given weights vector @binaryWeights
        """
        num_real_states = self.getNumPP()
        num_states = self.getNumStates()
        
        featuresList = self.binaryFeaturesList[:6]
        featuresMat = [[0 for j in range(num_real_states)] for i in range(num_real_states)]
        
        db = 0
        
        for sid1, state1 in enumerate(self.states[:num_real_states]):
            for sid2, state2 in enumerate(self.states[:num_real_states]):
                for i, f in enumerate(featuresList):
                    featureValue = f(sid1, state1, sid2, state2, db)
                    if featureValue == 1:
                        if featuresMat[sid1][sid2] == 0:
                            featuresMat[sid1][sid2] = 1
                        else:
                            raise Exception("WHAT")
                            
        print featuresMat
                        


    def viterbiPath(self, samples, M, P, MSC, PSC, mixture, inferHighLevelLabels=True):
        """
        Viterbi decoding of the most probable path.
        """
        num_states = self.getNumStates()
        num_real_states = self.getNumPP() 
        transitions = self.getEdgePotential(self.binaryWeights)
        distBin = self.distBin
        self.avgCoverage = self.estimateCoverage(samples)

        n = len(samples)
        predecessor = [[0 for i in range(num_states)] for j in xrange(n+1)] 
        #DP table dim: seq length+1 x num_states
        table = [[self.neg_inf for i in range(num_states)] for j in xrange(n+1)] 
        
        start_id, start_state = self.getStartState()
        exit_id, exit_state = self.getExitState()
        
        
        '''INITIALIZE'''
        #multiply the CNV prior into transition probabilities
        #self.adjustTransitionProbForPos(0, transitions)
        
        #the start state probability is 1 -> log(1)=0
        table[0][start_id] = 0. 
        #propagate over all silent states
        for state_id, state in enumerate(self.states[num_real_states:]):
            state_id += num_real_states
            if transitions[start_id][state_id][distBin[0]] == self.neg_inf: continue #just for speed-up
            table[0][state_id] = table[0][start_id] + transitions[start_id][state_id][distBin[0]]
            predecessor[0][state_id] = -1
         
        '''DYNAMIC PROGRAMMING'''
        #for all SNP positions do:
        for pos in xrange(1, n+1):
            #real positions are <1..n+1), pos 1 is the base case
            
            nodePotential = self.getNodePotential(pos, self.unaryWeights, samples, M, P, MSC, PSC, mixture)
            
            #(i) compute new values for all phased patterns - "real" states of the HMM:
            for state_id, state in enumerate(self.states[:num_real_states]):
                #emission probability in the given state
                #emis_p = self.logLHGivenStateWCoverage(pos-1, samples[pos-1], M[pos-1], P[pos-1], MSC[pos-1], PSC[pos-1], mixture, state)
                emis_p = nodePotential[state_id]
                
                #porobability of this state is `emission * max{previous state * transition}`
                max_prev = self.neg_inf
                arg_max = -1
                for prev_id in range(num_states):
                    if transitions[prev_id][state_id][distBin[pos]] == self.neg_inf: continue #just for speed-up
                    tmp = table[pos-1][prev_id] + transitions[prev_id][state_id][distBin[pos]]
                    if tmp > max_prev:
                        max_prev = tmp
                        arg_max = prev_id
                table[pos][state_id] = max_prev + emis_p
                predecessor[pos][state_id] = arg_max
            
            #multiply the CNV prior into transition probabilities
            #transitions = copy.deepcopy(self.transitions)
            #self.adjustTransitionProbForPos(pos, transitions)
            
            #(ii, iii) transitions from 'real' states and silent states with lower id to silent states
            #note: the states in self.states are already ordered such that a transistion from a silent 
            #   state to another silent state always points to a state with heigher index
            for state_id, state in enumerate(self.states[num_real_states:]):
                state_id += num_real_states
                #porobability of this state is `max{(actual real state or silent state with lower id) * transition}`
                max_prev = self.neg_inf
                arg_max = -1
                for prev_id in range(state_id + 1):
                    if transitions[prev_id][state_id][distBin[pos+1]] == self.neg_inf: continue #just for speed-up
                    tmp = table[pos][prev_id] + transitions[prev_id][state_id][distBin[pos+1]]
                    if tmp > max_prev:
                        max_prev = tmp
                        arg_max = prev_id
                table[pos][state_id] = max_prev
                predecessor[pos][state_id] = arg_max
            
            #normalize to sum to 1
            self.logNormalize(table[pos])
        
        '''RESULT'''
        path = [-1 for i in xrange(n+1)]
        path[n] = exit_id 
        while path[n] >= num_real_states: path[n] = predecessor[n][path[n]]
        for i in reversed(xrange(n)):
            path[i] = predecessor[i+1][path[i+1]]
            while path[i] >= num_real_states: path[i] = predecessor[i][path[i]]
        
        state_path = [[] for i in xrange(n+1)]
        for i in xrange(1, n+1):
            state = self.states[ path[i] ]
            if inferHighLevelLabels:
                path[i] = self.IPtoID[state.inheritance_pattern]
            state_path[i] = (state.inheritance_pattern, state.phased_pattern)
        
        print "Viterbi path probability: ", table[n][exit_id]    
        return path[1:], state_path[1:]
        

    def computeForward(self, samples, M, P, MSC, PSC, mixture):
        '''
        Posterior decoding: forward algorithm
        '''
        num_states = self.getNumStates()
        num_real_states = self.getNumPP()
        patterns = self.inheritance_patterns
        transitions = self.getEdgePotential(self.binaryWeights)
        distBin = self.distBin
        self.avgCoverage = self.estimateCoverage(samples)
        
        n = len(samples)
        #scale factors
        scale = [0. for i in xrange(n+1)]
        #DP table dim: seq length+1 x num_states
        table = [[self.neg_inf for i in range(num_states)] for j in xrange(n+1)] 
        
        start_id, start_state = self.getStartState()
        exit_id, exit_state = self.getExitState()
        
        
        '''INITIALIZE'''
        #multiply the CNV prior into transition probabilities
        #self.adjustTransitionProbForPos(0, transitions)
        
        #the start state probability is 1 -> log(1)=0
        table[0][start_id] = 0. 
#        #propagate over all silent states
#        for state_id, state in enumerate(self.states[num_real_states:]):
#            state_id += num_real_states
#            if transitions[start_id][state_id][distBin[0]] == self.neg_inf: continue #just for speed-up
#            table[0][state_id] = table[0][start_id] + transitions[start_id][state_id][distBin[0]]
        
        '''DYNAMIC PROGRAMMING'''
        #for all SNP positions do:
        for pos in xrange(1, n+1):
            #real positions are <1..n+1), pos 1 is the base case
            
            nodePotential = self.getNodePotential(pos, self.unaryWeights, samples, M, P, MSC, PSC, mixture)
            
            #(i) compute new values for all phased patterns - "real" states of the HMM:
            for state_id, state in enumerate(self.states[:num_real_states]):
                #emission probability in the given state
                #emis_p = self.logLHGivenStateWCoverage(pos-1, samples[pos-1], M[pos-1], P[pos-1], MSC[pos-1], PSC[pos-1], mixture, state)
                emis_p = nodePotential[state_id]
                
                #probability of this state is `emission * \sum {previous state * transition}`
                summ = self.neg_inf
                for prev_id in range(num_states):
                    if transitions[prev_id][state_id][distBin[pos]] == self.neg_inf: continue #just for speed-up
                    tmp =  table[pos-1][prev_id] + transitions[prev_id][state_id][distBin[pos]]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = emis_p + summ
            
            #multiply the CNV prior into transition probabilities
            #transitions = copy.deepcopy(self.transitions)
            #self.adjustTransitionProbForPos(pos, transitions)
            
            #(ii, iii) transitions from 'real' states and silent states with *lower* id to silent states
            #note: the states in self.states are already ordered such that a transistion from a silent 
            #   state to another silent state always points to a state with heigher index
            for state_id, state in enumerate(self.states[num_real_states:]):
                state_id += num_real_states
                
                #probability of this state is `\sum {(actual real state or silent state with lower id) * transition}`
                summ = self.neg_inf
                for prev_id in range(state_id):
                    if transitions[prev_id][state_id][distBin[pos+1]] == self.neg_inf: continue #just for speed-up
                    tmp =  table[pos][prev_id] + transitions[prev_id][state_id][distBin[pos+1]]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = summ
                
            #normalize to sum to 1
            scale[pos] = self.logNormalize(table[pos])
        
        '''RESULT'''
        #p(X) - probability of the observed sequence (samples)
        pX = table[n][exit_id]
        #print "fwd: ", math.exp(pX), pX

        return table, pX, scale
        
        
    def computeBackward(self, samples, M, P, MSC, PSC, mixture, scale):
        '''
        Posterior decoding: backward algorithm
        '''
        num_states = self.getNumStates()
        num_real_states = self.getNumPP()
        patterns = self.inheritance_patterns
        transitions = self.getEdgePotential(self.binaryWeights)
        distBin = self.distBin
        self.avgCoverage = self.estimateCoverage(samples)
        
        n = len(samples)
        #DP table dim: seq length+2 x num_states
        table = [[self.neg_inf for i in range(num_states)] for j in xrange(n+1)]
        
        start_id, start_state = self.getStartState()
        exit_id, exit_state = self.getExitState()
        
        
        '''INITIALIZE'''
        #multiply the CNV prior into transition probabilities
        #self.adjustTransitionProbForPos(n-1, transitions)
        
        #the exit state probability is 1 -> log(1)=0
        table[n][exit_id] = 0. 
        #(iii) add transitions from silent to silent states with *higher* id (*no emission*)
        for state_id in reversed(range(num_real_states, num_states)):
            #\sum {'next'_state * transition_from_here} (no emission)
            summ = self.neg_inf
            for next_id in reversed(range(state_id, num_states)):
                if transitions[state_id][next_id][distBin[n+1]] == self.neg_inf: continue #just for speed-up
                tmp = transitions[state_id][next_id][distBin[n+1]] + table[n][next_id]
                summ = self.logSum(summ, tmp)
            table[n][state_id] = self.logSum(summ, table[n][state_id])
        
        #(ii) add transitions from 'real' states to silent states (*no emission*)
        for state_id in range(num_real_states):
            #\sum {'next'_state * transition_from_here} (no emission)
            summ = self.neg_inf
            for next_id in range(num_real_states, num_states):
                if transitions[state_id][next_id][distBin[n+1]] == self.neg_inf: continue #just for speed-up
                tmp = transitions[state_id][next_id][distBin[n+1]] + table[n][next_id]
                summ = self.logSum(summ, tmp)
            table[n][state_id] = self.logSum(summ, table[n][state_id])
        
        '''DYNAMIC PROGRAMMING'''    
        #for all SNP positions do:
        for pos in reversed(xrange(0, n)):
            #real positions are <1..n+1), the 0 and n+1 are sentinels
            
            #precompute emission probabilities
            #emis_p = []
            #for state in self.states[:num_real_states]:
            #    emis_p.append(self.logLHGivenStateWCoverage(pos, samples[pos], M[pos], P[pos], MSC[pos], PSC[pos], mixture, state))
            nodePotential = self.getNodePotential(pos+1, self.unaryWeights, samples, M, P, MSC, PSC, mixture)
            
            #(i) transitions from all states to 'real' states of the HMM:
            for state_id in range(num_states):
                # \sum {emission_in_'next' * 'next'_state * transition_from_here}
                summ = self.neg_inf
                for next_id in range(num_real_states):
                    if transitions[state_id][next_id][distBin[pos+1]] == self.neg_inf: continue #just for speed-up
                    tmp = transitions[state_id][next_id][distBin[pos+1]] + table[pos+1][next_id] + nodePotential[next_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = summ
            
            #(iii) add transitions from silent to silent states with *higher* id (*no emission*)
            for state_id in reversed(range(num_real_states, num_states)):
                #\sum {'next'_state * transition_from_here} (no emission)
                summ = self.neg_inf
                for next_id in reversed(range(state_id, num_states)):
                    if transitions[state_id][next_id][distBin[pos+1]] == self.neg_inf: continue #just for speed-up
                    tmp = transitions[state_id][next_id][distBin[pos+1]] + table[pos][next_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = self.logSum(summ, table[pos][state_id])
            
            #(ii) add transitions from 'real' states to silent states (*no emission*)
            for state_id in range(num_real_states):
                #\sum {'next'_state * transition_from_here} (no emission)
                summ = self.neg_inf
                for next_id in range(num_real_states, num_states):
                    if transitions[state_id][next_id][distBin[pos+1]] == self.neg_inf: continue #just for speed-up
                    tmp = transitions[state_id][next_id][distBin[pos+1]] + table[pos][next_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = self.logSum(summ, table[pos][state_id])
            
            #multiply the CNV prior into transition probabilities
            #transitions = copy.deepcopy(self.transitions)
            #self.adjustTransitionProbForPos(pos-1, transitions)
            
            #pseudonormalize by scaling factor from fwd
            for state_id in range(num_states):
                table[pos][state_id] -= scale[pos+1] 
        
        '''RESULT'''
        #p(X) - probability of the observed sequence (samples)
        pX = table[0][start_id]
        #print "bck: ", math.exp(pX), pX
        
        return table, pX
        
        
    def maxPosteriorDecoding(self, samples, M, P, MSC, PSC, mixture):
        """Maximum posterior probability decoding. 
        
        Returns list of states -- for each position the one with highest posterior
        probability.
        
        """
        n = len(samples)
        
        fwd, pX1, scale = self.computeForward(samples, M, P, MSC, PSC, mixture)
        bck, pX2 = self.computeBackward(samples, M, P, MSC, PSC, mixture, scale)
        
        if abs(pX1-pX2) > 10e-8:
            print "p(X) from forward and backward DP doesn't match", pX1, pX2
            return
        
        path = []
        for pos in xrange(1, n+1):
            max_posterior = self.neg_inf
            arg_max = -1
            for s in range(self.getNumPP()):
                # fwd, bck are in log space (therefore + instead *); p(X) is constant
                tmp = fwd[pos][s] + bck[pos][s]
                if tmp > max_posterior:
                    max_posterior = tmp
                    arg_max = s
            #print max_posterior
            max_IP = self.states[ arg_max ].inheritance_pattern
            path.append( self.IPtoID[max_IP] )
        
        return path
    
    
    def posteriorDecoding(self, samples, M, P, MSC, PSC, mixture):
        """Posterior probability decoding.
        
        For each position returns a list of states ordered by their posterior
        porobability. Returns list(list(tuple(probability, state))) and log of p(X)
        
        """
        n = len(samples)
        num_real_states = self.getNumPP()
        
        fwd, pX1, scale = self.computeForward(samples, M, P, MSC, PSC, mixture)
        bck, pX2 = self.computeBackward(samples, M, P, MSC, PSC, mixture, scale)
        
        logZ = pX1 + sum(scale)
        #print pX1, pX2, sum(scale), Z
        
        if abs(pX1-pX2) > 10e-8:
            print "p(X) from forward and backward DP doesn't match", pX1, pX2
            return
        
        table = []
        over90 = 0
        numall = 0
        for pos in xrange(1, n+1):
            tmp = []
            tmp2 = []
            for ip_id, ip in enumerate(self.inheritance_patterns):
                #max over corresponding phased patterns
                max_ = self.neg_inf
                for s_id, s in enumerate(self.states[:num_real_states]):
                    if s.inheritance_pattern == ip:
                        #fwd, bck are in log space (therefore + instead *); p(X) is constant
                        #max_ = max(max_, fwd[pos][s_id] + bck[pos][s_id])
                        max_ = max(max_, fwd[pos][s_id] + bck[pos][s_id] - pX1)
                
                '''#sum over corresponding phased patterns
                sum_ = self.neg_inf
                for s_id, s in enumerate(self.states):
                    if s.inheritance_pattern == ip:
                        #fwd, bck are in log space (therefore + instead *); p(X) is constant
                        sum_ = self.logSum(sum_, fwd[pos][s_id] + bck[pos][s_id])
                max_ = sum_
                '''
                    
                tmp.append((max_, ip_id))
                tmp2.append([max_, ip_id])
            
            sum_ = self.neg_inf
            for x in tmp2:
                sum_ = self.logSum(sum_, x[0])
            for i in range(len(tmp2)):
                tmp2[i][0] -= sum_
                tmp2[i][0] = int(math.exp(tmp2[i][0])*1000)
            tmp2.sort(reverse=True)
            #print tmp2
            if tmp2[0][0] >= 900: over90 += 1
            numall += 1
            
            tmp.sort(reverse=True)
            table.append(tmp)
        
        print over90, numall, "that is:", over90/float(numall)
        return table, logZ
        
            
    def likelihoodDecoding(self, samples, M, P, MSC, PSC, mixture):
        '''
        for each position returns a list of states ordered by the likelihood
        of the data in corresponding state, i.e. by emission porobabilities
        '''
        n = len(samples)
        num_real_states = self.getNumPP()
        
        table = []
        for pos in xrange(1, n+1):
            tmp = []
            for ip_id, ip in enumerate(self.inheritance_patterns):
                #max over corresponding phased patterns
                max_ = self.neg_inf
                for s in self.states[:num_real_states]:
                    if s.inheritance_pattern == ip:
                        ll = self.logLHGivenStateWCoverage(pos-1, samples[pos-1], M[pos-1], P[pos-1], MSC[pos-1], PSC[pos-1], mixture, s)
                        max_ = max(max_, ll)
                
                '''#sum over corresponding phased patterns
                sum_ = self.neg_inf
                for s in self.states:
                    if s.inheritance_pattern == ip:
                        ll = self.logLHGivenState(samples[pos-1], M[pos-1], P[pos-1], MSC[pos-1], PSC[pos-1], mixture, s)
                        sum_ = self.logSum(sum_, ll)
                '''        

                tmp.append((max_ , ip_id))
                
            tmp.sort(reverse=True)    
            table.append(tmp)
        return table
                        
    
    def maxLikelyState(self, nuc_counts, M, P, MSC, PSC, mixture):
        '''
        blah blah
        import fcnvHMM; f = fcnvHMM.FCNV()
        tmp = f.maxLikelyState((100,0,0,100), ('A','T'), ('A','T'), 0.1)

        '''
        num_real_states = self.getNumPP()
        
        tmp = []
        for ip_id, ip in enumerate(self.inheritance_patterns):
            #max over corresponding phased patterns
            max_ = self.neg_inf
            for s in self.states[:num_real_states]:
                if s.inheritance_pattern == ip:
                    ll = self.logLHGivenState(nuc_counts, M, P, MSC, PSC, mixture, s)
                    ll2 =self.logLHGivenState(nuc_counts, M, P, MSC, PSC, mixture, s, True)
                    max_ = max(max_, ll)

                    tmp.append((ll, ll2, s.phased_pattern, s.inheritance_pattern))
            
        tmp.sort(reverse=True)    

        return tmp

if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    
############# edges when using "silent states" 
#def getEdgePotential(self, binaryWeights):
#        """
#        Compute overall edge potential shared across all edges
#        
#        return 'num states' x 'num states' matrix representing the potential that
#        pools all edge features together given weights vector @binaryWeights
#        """
#        num_real_states = self.getNumPP()
#        num_states = self.getNumStates()
#        
#        edgePot = [[0. for x in range(num_states)] for y in range(num_states)]
#        
#        wStayNormal = binaryWeights[0]
#        wRecombNormal = binaryWeights[1]
#        wGoNormal = binaryWeights[2]
#        wStayCNV = binaryWeights[3]
#        wRecombCNV = binaryWeights[4]
#        wGoCNV = binaryWeights[5]
#        wEps = binaryWeights[6]
#        
#        #first generate pairwise energy from the real states
#        for ip_id, ip in enumerate(self.inheritance_patterns):
#            
#            if ip == (1, 1): #generate Normal states pairwise energy
#                for i, state1 in enumerate(self.states[:num_real_states]):
#                    if state1.inheritance_pattern != ip: continue
#                    #states "inside the IP component"
#                    for j, state2 in enumerate(self.states[:num_real_states]):
#                        if i == j: #stay in the state
#                            edgePot[i][i] = wStayNormal
#                        elif state2.inheritance_pattern == ip: #to recombination of the same IP
#                            edgePot[i][j] = wRecombNormal
#                    #to the silent exit node    
#                    outState_id= self.getCorrespondingOutState(state1)[0]
#                    edgePot[i][outState_id] = wGoNormal
#                    
#            else: #generate CNV states pairwise energy
#                for i, state1 in enumerate(self.states[:num_real_states]):
#                    if state1.inheritance_pattern != ip: continue
#                    #states "inside the IP component"
#                    for j, state2 in enumerate(self.states[:num_real_states]):
#                        if i == j: #stay in the state
#                            edgePot[i][i] = wStayCNV
#                        elif state2.inheritance_pattern == ip: #to recombination of the same IP
#                            if state1.inheritance_pattern == (2, 1):
#                                if state1.phased_pattern[0][0] ^ state1.phased_pattern[0][1] != state2.phased_pattern[0][0] ^ state2.phased_pattern[0][1]:
#                                    edgePot[i][j] = wEps
#                                    continue
#                            if state1.inheritance_pattern == (1, 2):
#                                if state1.phased_pattern[1][0] ^ state1.phased_pattern[1][1] != state2.phased_pattern[1][0] ^ state2.phased_pattern[1][1]:
#                                    edgePot[i][j] = wEps
#                                    continue
#                            edgePot[i][j] = wRecombCNV
#                    
#                    #to the silent exit node    
#                    outState_id = self.getCorrespondingOutState(state1)[0]
#                    edgePot[i][outState_id] = wGoCNV
#                    
#        #now generate pairwise energy from the silent states
#        inNormal_id = self.getCorrespondingInState( HMMState((1,1), ()) )[0]
#        constW = 0.2
#        for i, state1 in enumerate(self.states[num_real_states:]):
#            i += num_real_states
#            #if it is the start node
#            if state1.phased_pattern == "s":
#                for j, state2 in enumerate(self.states[num_real_states:]):
#                    if state2.phased_pattern == "in":
#                        edgePot[i][num_real_states + j] = constW
#            
#            #if it is a silent exit node
#            elif state1.phased_pattern == "out":
#                j = self.getExitState()[0]
#                edgePot[i][j] = constW
#                        
#                if self.isNormal(state1):
#                    for j, state2 in enumerate(self.states[num_real_states:]):
#                        if state2.phased_pattern == "in" and \
#                         state2.inheritance_pattern != state1.inheritance_pattern \
#                         and state2.inheritance_pattern != (0,2) \
#                         and state2.inheritance_pattern != (2,0):
#                            edgePot[i][num_real_states + j] = constW
#                else:
#                    edgePot[i][inNormal_id] = constW
#                    
#            #if it is a silent starting node
#            elif state1.phased_pattern == "in":
#                for j, state2 in enumerate(self.states[:num_real_states]):
#                    if state2.inheritance_pattern == state1.inheritance_pattern:
#                        edgePot[i][j] = constW      
#        
#        #take logarithm
#        for i in range(num_states):
#            for j in range(num_states):
#                if edgePot[i][j] < 10e-10: edgePot[i][j] = self.neg_inf
#                else: edgePot[i][j] = math.log(edgePot[i][j])
#        
#        return edgePot
    
    
