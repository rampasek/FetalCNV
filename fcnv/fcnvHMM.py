import random
import math
import numpypy as np
import itertools
import copy


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
    
    def __init__(self):
        """Initialize new FCNV object"""
        super(FCNV, self).__init__()
        
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
        
        #generate inheritance patterns
        self.inheritance_patterns = []
        for mats in range(3):
            for pats in range(3):
                if mats+pats in [0,4]: continue
                self.inheritance_patterns.append((mats, pats))
        
               
        #generate phased variants of the inheritance patterns and list of HMM States
        self.states = []
        self.phased_patterns = []
        for ip in self.inheritance_patterns:
            for mPhased in itertools.combinations_with_replacement([0,1], ip[0]):
                for pPhased in itertools.combinations_with_replacement([0,1], ip[1]):
                    phased_pattern = (tuple(mPhased), tuple(pPhased))
                    self.phased_patterns.append(phased_pattern)
                    new_state = HMMState(ip, phased_pattern)
                    self.states.append(new_state)
        
        self.main_state = 3
        self.main_states = []
        for i, s in enumerate(self.states):
            if s.inheritance_pattern == (1,1): self.main_states.append(i)
        
        #generate the table of transition probaboilities
        self.transitions = self.generateTransitionProb()

        #print the table
        num_states = self.getNumStates()
        for k in xrange(num_states):
            for l in xrange(num_states):
                print "%0.4f" % (math.exp(self.transitions[k][l])),
            print " "
            
    
    def isNormal(self, state):
        return state.inheritance_pattern == (1,1)
    
    def areRecombination(self, state1, state2):
        if state1.inheritance_pattern != state2.inheritance_pattern:
            return False
        eq = 0
        for hap in range(2):
            eq += (state1.phased_pattern[hap] == state2.phased_pattern[hap])
        return (eq > 0)
    
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
    
    def generateTransitionProb(self):
        """Generate transition probabilities between HMM states"""
        #number of possible states per one SNP position
        num_states = self.getNumStates()
        
        #first generate pseudocounts
        trans = [[0. for x in range(num_states)] for y in range(num_states)]
        for i, s in enumerate(self.states):
            for j, t in enumerate(self.states):
                if s == t: #remain in this state
                    trans[i][j] = 0.999
                    if self.isNormal(s): trans[i][j] = 0.9999
                elif self.isNormal(s) and self.isNormal(t): #recombination in normal IP
                    trans[i][j] += 0.0001
                elif self.isNormal(s) ^ self.isNormal(t): #CNV <-> normal IP
                    if self.areAdjacent(s, t) \
                       or s.inheritance_pattern in [(2,0), (0,2)] \
                       or t.inheritance_pattern in [(2,0), (0,2)]:
                        trans[i][j] += 0.001
                elif self.areRecombination(s, t): #recombination within CNV
                    trans[i][j] += 0.00001
                    
        
        #normalize and take logarithm
        for i in range(num_states):
            for j in range(num_states):
                if trans[i][j] < 10e-8: trans[i][j] = float("-inf")
                else: trans[i][j] = math.log(trans[i][j])
            self.logNormalize(trans[i])

        return trans
    
    def getNumIP(self):
        """Return number of inheritance patterns"""
        return len(self.inheritance_patterns)
        
    def getNumPP(self):
        """Return number of phased inheritance patterns"""
        return len(self.phased_patterns)
        
    def getNumStates(self):
        """Return number of HMM states per one SNP position"""
        return self.getNumPP()
    
    def logNormalize(self, l):
        """Normalize the given list of log-probabilities.
        
        Normalizes the list in place and returns the used scale factor.
        
        """
        sum_ = "nothing"
        for x in l: sum_ = self.logSum(sum_, x)
        for i in range(len(l)): l[i] -= sum_
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
        
        if x == "nothing": return y
        elif y == "nothing": return x
        elif x == float("-inf") and y == float("-inf"): return float("-inf")
        elif x == float("inf") and y == float("inf"): return float("inf")
        
        a = max(x, y)
        b = min(x, y)
        if a == float("inf") and b == float("-inf"): return float("nan")
        
        #return a + lookup(b-a):
        return a + math.log(1 + math.exp(b-a) )
    
    def logMultinomial(self, xs, ps):
        """Calculate probability mass function of multinomial distribution
        returns log(Multi(xs | sum(xs), ps))
        
        Arguments:
        xs -- [ints] -- counts
        ps -- [floats] -- probabilities
        returns float 
        
        >>> f = FCNV()
        >>> f.logMultinomial([1,2,3], [.2,.3,.5])
        -2.0024805005437072
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
            if isinstance(x, list):
                res = []
                for val in x:
                    res.append(gammaln(val+1))
                return res
            else: 
                return gammaln(x+1)
        
        n = sum(xs)
        '''#numpy implementation:
        xs, ps = np.array(xs), np.array(ps)
        result = logFactorial(n) - sum(logFactorial(xs)) + sum(xs * np.log(ps))
        '''
        
        result = logFactorial(n) - sum(logFactorial(xs))
        for i in range(len(ps)):
            result += xs[i] * math.log(ps[i])
            
        return result

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
        if val: return val
        '''

        adjusted_fetal_admix = mix/2. * len(fetal)
        adjusted_maternal_admix = (1.-mix) / 2. * len(maternal)
        cmix = adjusted_fetal_admix / (adjusted_maternal_admix + adjusted_fetal_admix)
        dist = []
        for nuc in self.nucleotides:
            p = maternal.count(nuc) / float(len(maternal)) * (1.-cmix)
            p += fetal.count(nuc) / float(len(fetal)) * (cmix)
            dist.append(p)
        
        summ = 0.
        for x in dist:
            summ += 0.01 + x #add noise
        for i in range(len(dist)):
            dist[i] = (dist[i] + 0.01) / summ #normalize
        
        '''self.distributionCache[code] = dist'''
        return dist
    
    
    #TODO: NEEDS OPTIMIZATION if possible !!!!!!!!!!!!!!!!!    
    def logLHGivenState(self, nuc_counts, maternal_alleles, paternal_alleles, mix, state):
        '''
        >>> f = FCNV()
        >>> f.logLHGivenState([3, 0, 0, 71], ['T', 'T'], ['T', 'A'], 0.1, [1,1])
        -3.6974187244239025
        >>> f.logLHGivenState([3, 0, 0, 71], ['T', 'T'], ['T', 'A'], 0.1, [2,1]) #-3.4838423152150568
        -3.6507523341244412
        >>> f.logLHGivenState([3, 0, 0, 71], ['T', 'T'], ['T', 'A'], 0.1, [0,2])
        -3.1614953504205028
        '''
        pattern = state.phased_pattern
               
        fetal_alleles = []
        for mHpatt in pattern[0]:
            fetal_alleles.append(maternal_alleles[mHpatt])
        for pHpatt in pattern[1]:
            fetal_alleles.append(paternal_alleles[pHpatt])
        fetal_alleles.sort()
        
        #Caching:
        code = (tuple(nuc_counts), tuple(maternal_alleles), tuple(fetal_alleles), mix)
        val = self.logLikelihoodCache.get(code)
        if val: return val
        
        ps = self.allelesDistribution(maternal_alleles, fetal_alleles, mix)
        result = self.logMultinomial(nuc_counts, ps)
        
        self.logLikelihoodCache[code] = result
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
        for i in xrange(len(samples)):
            if M[i][0] == M[i][1] and not P[i][0] == P[i][1] and \
              (P[i][0] == M[i][0] or P[i][1] == M[i][0]):
                #print M[i], P[i], samples[i]            
                type1 = M[i][0]
                type2 = P[i][0]
                if type1 == type2: type2 = P[i][1]
                num_type2 = samples[i][self.nucleotides.index(type2)]

                if not num_type2 == 0:
                    sum_all = sum(samples[i])
                    num += 1
                    mix += (2.*num_type2)/sum_all
        return mix/num
    
    
    def viterbiPath(self, samples, M, P, mixture):
        """Viterbi decoding of the most probable path.
        
        """
        num_states = self.getNumStates()
        transitions = self.transitions
        
        n = len(samples)
        predecessor = [[0 for i in range(num_states)] for j in xrange(n+1)] 
        #DP table dim: seq length+1 x num_states
        table = [[0. for i in range(num_states)] for j in xrange(n+1)] 
        
        #initital base case
        for state_id, state in enumerate(self.states):
            emis_p = self.logLHGivenState(samples[0], M[0], P[0], mixture, state)
            table[1][state_id] = math.log(1./num_states) + emis_p #initial state probabilities   
            predecessor[1][state_id] = -1
        
        #for all SNP positions do:
        for pos in xrange(2, n+1):
            #real positions are <1..n+1), pos 1 is the base case

            #for each SNP consider all inheritance patterns - states of the HMM:
            for state_id, state in enumerate(self.states):
                #emission probability in the given state
                emis_p = self.logLHGivenState(samples[pos-1], M[pos-1], P[pos-1], mixture, state)
                
                #porobability of this state is emission * max{previous state * transition}
                max_prev = float("-inf")
                arg_max = -1
                for prev_id in range(num_states):
                    if transitions[prev_id][state_id] == float("-inf"): continue #just for speed-up
                    tmp = table[pos-1][prev_id] + transitions[prev_id][state_id]
                    if tmp > max_prev:
                        max_prev = tmp
                        arg_max = prev_id
                table[pos][state_id] = max_prev + emis_p
                predecessor[pos][state_id] = arg_max
            
            #normalize to sum to 1
            self.logNormalize(table[pos])
        
        path = [-1 for i in xrange(n+1)]
        maxx = max(table[n])
        path[n] = table[n].index(maxx)
        for i in reversed(xrange(n)):
            path[i] = predecessor[i+1][path[i+1]]
        
        for i in xrange(1, n+1):
            state = self.states[ path[i] ]
            path[i] = self.inheritance_patterns.index(state.inheritance_pattern)
            
        return path[1:]
    
    def computeForward(self, samples, M, P, mixture):
        '''
        Posterior decoding: forward algorithm
        '''
        num_states = self.getNumStates()
        patterns = self.inheritance_patterns
        transitions = self.transitions
        
        n = len(samples)
        #scale factors
        scale = [0. for i in xrange(n+1)]
        #DP table dim: seq length+1 x num_states
        table = [[float("-inf") for i in range(num_states)] for j in xrange(n+1)] 

        #initital base case
        for state_id, state in enumerate(self.states):
            emis_p = self.logLHGivenState(samples[0], M[0], P[0], mixture, state)
            table[1][state_id] = math.log(1./num_states) + emis_p #initial state probabilities

        #for all SNP positions do:
        for pos in xrange(2, n+1):
            #real positions are <1..n+1), pos 1 is the base case

            #for each SNP consider all inheritance patterns - states of the HMM:
            for state_id, state in enumerate(self.states):
                #emission probability in the given state
                emis_p = self.logLHGivenState(samples[pos-1], M[pos-1], P[pos-1], mixture, state)
                
                #probability of this state is emission * \sum {previous state * transition}
                summ = "nothing"
                for prev_id in range(num_states):
                    if transitions[prev_id][state_id] == float("-inf"): continue #just for speed-up
                    tmp =  table[pos-1][prev_id] + transitions[prev_id][state_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = emis_p + summ
            
            #normalize to sum to 1
            scale[pos] = self.logNormalize(table[pos])
            
        #p(X) - probability of the observed sequence (samples)
        pX = "nothing"
        for i in table[n]:
            pX = self.logSum(pX, i)
        #print "fwd: ", math.exp(pX), pX

        return table, pX, scale
        
    def computeBackward(self, samples, M, P, mixture, scale):
        '''
        Posterior decoding: backward algorithm
        '''
        num_states = self.getNumStates()
        patterns = self.inheritance_patterns
        transitions = self.transitions
        
        n = len(samples)
        #DP table dim: seq length+2 x num_states
        table = [[float("-inf") for i in range(num_states)] for j in xrange(n+1)]
        
        #initital base case
        for state_id in range(num_states):
            table[n][state_id] = 0. #initial state probabilities are 1 -> log(1)=0
            
        #for all SNP positions do:
        for pos in reversed(xrange(1, n)):
            #real positions are <1..n+1), the 0 and n+1 are sentinels
            
            #precompute emission probabilities
            emis_p = []
            for state in self.states:
                emis_p.append(self.logLHGivenState(samples[pos], M[pos], P[pos], mixture, state) )
            
            #for each SNP consider all inheritance patterns - states of the HMM:
            for state_id in range(num_states):
                #probability of this state is:
                # \sum {emission_in_'next' * 'next'_state * transition_from_here}
                summ = "nothing"
                for next_id in range(num_states):
                    if transitions[state_id][next_id] == float("-inf"): continue #just for speed-up
                    tmp = transitions[state_id][next_id] + table[pos+1][next_id] + emis_p[next_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = summ - scale[pos+1] #pseudonormalize by scaling factor from fwd
            
            
        #p(X) - probability of the observed sequence (samples)
        pX = "nothing"
        for s in range(num_states):
            emis_p = self.logLHGivenState(samples[0], M[0], P[0], mixture, self.states[s])
            tmp = math.log(1./num_states) + table[1][s] + emis_p - scale[1] #scale by factor from fwd
            pX = self.logSum(pX, tmp)
        #print "bck: ", math.exp(pX), pX
        
        return table, pX
        
        
    def maxPosteriorDecoding(self, samples, M, P, mixture):
        """Maximum posterior probability decoding. 
        
        Returns list of states -- for each position the one with highest posterior
        probability.
        
        """
        n = len(samples)
        
        fwd, pX1, scale = self.computeForward(samples, M, P, mixture)
        bck, pX2 = self.computeBackward(samples, M, P, mixture, scale)
        
        if abs(pX1-pX2) > 10e-8:
            print "p(X) from forward and backward DP doesn't match", pX1, pX2
            return
        
        path = []
        for pos in xrange(1, n+1):
            max_posterior = float("-inf")
            arg_max = -1
            for s in range(self.getNumStates()):
                # fwd, bck are in log space (therefore + instead *); p(X) is constant
                tmp = fwd[pos][s] + bck[pos][s]
                if tmp > max_posterior:
                    max_posterior = tmp
                    arg_max = s
            #print max_posterior
            max_IP = self.states[ arg_max ].inheritance_pattern
            path.append( self.inheritance_patterns.index(max_IP) )
        
        return path
    
    
    def posteriorDecoding(self, samples, M, P, mixture):
        """Posterior probability decoding.
        
        For each position returns a list of states ordered by their posterior
        porobability. Returns list(list(tuple(probability, state)))
        
        """
        n = len(samples)
        
        fwd, pX1, scale = self.computeForward(samples, M, P, mixture)
        bck, pX2 = self.computeBackward(samples, M, P, mixture, scale)
        
        if abs(pX1-pX2) > 10e-8:
            print "p(X) from forward and backward DP doesn't match", pX1, pX2
            return
        
        table = []
        for pos in xrange(1, n+1):
            tmp = []
            for ip_id, ip in enumerate(self.inheritance_patterns):
                #max over corresponding phased patterns
                max_ = float("-inf")
                for s_id, s in enumerate(self.states):
                    if s.inheritance_pattern == ip:
                        #fwd, bck are in log space (therefore + instead *); p(X) is constant
                        max_ = max(max_, fwd[pos][s_id] + bck[pos][s_id])
                        
                tmp.append((max_, ip_id))
                
            tmp.sort(reverse=True)    
            table.append(tmp)
            
        return table
        
        
    #TODO: NEEDS REWRITING !!!!!!!!!!!!!!!!!    
    def baumWelshForTransitions(self, samples, M, P, mixture):
        """Baum-Welsh EM algorithm for HMM parameter estimation.
        
        Training implemented only for transition probabilities.
        
        Doesn't return anything, directly modifies self.transitions
        
        """
        num_states = self.getNumStates()
        n = len(samples)
        A = [[float("-inf") for x in range(num_states)] for y in range(num_states)]
        num_iter = 0
        while True and num_iter < 5:
            num_iter += 1
            #expectation
            fwd, pX1, scale = self.computeForward(samples, M, P, mixture)
            bck, pX2 = self.computeBackward(samples, M, P, mixture, scale)
            if abs(pX1-pX2) > 10e-8:
                print "p(X) from forward and backward DP doesn't match", pX1, pX2
                return  
            
            
            for k in range(num_states):            
                for l in range(num_states):
                    if self.transitions[k][l] == float("-inf"):
                        A[k][l] = float("-inf")
                        continue
                        
                    val = "nothing"
                    for pos in xrange(1, n):
                        emis_p = self.logLHGivenState(samples[pos], M[pos],\
                         P[pos], mixture, self.states[l])
                        tmp = fwd[pos][k] + self.transitions[k][l] + emis_p + bck[pos+1][l]
                        if tmp == float("-inf"): continue
                        val = self.logSum(val, tmp)
                    
                    if val == "nothing":
                        A[k][l] = float("-inf")
                    else:
                        val -= pX1
                        A[k][l] = val
                    #if math.exp(A[k][l]) < 0.1:
                    #    A[k][l] = math.log(0.1)
                         
            #maximization
            old_trans = copy.deepcopy(self.transitions)
            for k in range(num_states):
                summ = "nothing"
                for m in range(num_states):
                    if A[k][m] == float("-inf"): continue
                    summ = self.logSum(summ, A[k][m])
                #print summ, math.exp(summ)
                for l in range(num_states):
                    self.transitions[k][l] = A[k][l] - summ
                    #if self.transitions[k][l] < math.log(0.0005):
                    #    self.transitions[k][l] = math.log(0.0005)
            
            change = 0.
            for k in range(num_states):
                summ = "nothing"
                for m in range(num_states):
                    if A[k][m] == float("-inf"): continue
                    summ = self.logSum(summ, self.transitions[k][m])
                #print summ
                for l in range(num_states):
                    #self.transitions[k][l] -= summ
                    if not (self.transitions[k][l] == float("-inf") or old_trans[k][l] == float("-inf")):
                        change += (self.transitions[k][l] - old_trans[k][l])**2
                    
            print "change: ", change 
            if change < 1: break
                    
        #print the table    
        for k in xrange(num_states):
            for l in xrange(num_states):
                print "%0.4f" % (math.exp(self.transitions[k][l])), 
            print " "
    
    
    #TODO: NEEDS REWRITING !!!!!!!!!!!!!!!!!         
    def viterbiTrainingForTransitions(self, samples, M, P, mixture):
        """Viterbi training algorithm for HMM parameter estimation.
        
        Training implemented only for transition probabilities.
        
        Doesn't return anything, directly modifies self.transitions
        
        """
        num_states = self.getNumStates()
        n = len(samples)
        A = [[0. for x in range(num_states)] for y in range(num_states)]
        num_iter = 0
        while True and num_iter < 5:
            num_iter += 1
            #expectation
            path, table = self.viterbiPath(samples, M, P, mixture)
            p_path = max(table[len(samples)])
            print "P of the max probable path: ", p_path
            
            for k in range(num_states):
                for l in range(num_states):
                    val = 0
                    for pos in xrange(1, n-1):
                        if path[pos] == k and path[pos+1] == l:
                            val += 1
                    #if val < 1: val = 1
                    A[k][l] = val + int(k == self.main_state)
                A[k][k] += 1
                A[k][self.main_state] += 1 
                
                    
                         
            #maximization
            old_trans = copy.deepcopy(self.transitions)
            for k in range(num_states):
                summ = sum(A[k])
                if summ == 0:
                    summ = float("inf")
                else: 
                    summ = math.log(sum(A[k]))
                #for m in range(num_states):
                #    summ = self.logSum(summ, A[k][m])
                #print summ, math.exp(summ)
                for l in range(num_states):
                    if A[k][l] == 0:
                        self.transitions[k][l] = float("-inf")
                    else:
                        self.transitions[k][l] = math.log(A[k][l]) - summ
                    #if self.transitions[k][l] < math.log(0.0005):
                    #    self.transitions[k][l] = math.log(0.0005)
                        
            change = 0.
            for k in range(num_states):
                summ = "nothing"
                for m in range(num_states):
                    if A[k][m] == float("-inf"): continue
                    summ = self.logSum(summ, self.transitions[k][m])
                for l in range(num_states):
                    #self.transitions[k][l] -= summ
                    if not (self.transitions[k][l] == float("-inf") or old_trans[k][l] == float("-inf")):
                        change += (self.transitions[k][l] - old_trans[k][l])**2
                    
            print "change: ", change 
            if change < 1: break
                   
        #print the table    
        for k in xrange(num_states):
            for l in xrange(num_states):
                print "%0.4f" % (math.exp(self.transitions[k][l])), 
            print " "
        
            
    def likelihoodDecoding(self, samples, M, P, mixture):
        '''
        for each position returns a list of states ordered by the likelihood
        of the data in corresponding state, i.e. by emission porobabilities
        '''
        n = len(samples)
        
        table = []
        for pos in xrange(1, n+1):
            tmp = []
            for ip_id, ip in enumerate(self.inheritance_patterns):
                #max over corresponding phased patterns
                max_ = float("-inf")
                for s in self.states:
                    if s.inheritance_pattern == ip:
                        ll = self.logLHGivenState(samples[pos-1], M[pos-1], P[pos-1], mixture, s)
                        max_ = max(max_, ll)
                
                tmp.append((max_ , ip_id))
                
            tmp.sort(reverse=True)    
            table.append(tmp)
        return table

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
