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
        
        #generate silent states
        self.states.append( HMMState("s", "s") ) #start state
        
        normalIP = (1,1)
        for ip in self.inheritance_patterns:
            if ip == normalIP: continue
            new_state = HMMState(ip, "out")
            self.states.append(new_state)
        
        self.states.append( HMMState(normalIP, "in") )
        self.states.append( HMMState(normalIP, "out") )
        
        for ip in self.inheritance_patterns:
            if ip == normalIP: continue
            new_state = HMMState(ip, "in")
            self.states.append(new_state)
        
        self.states.append( HMMState("t", "t") ) #end state
        
        #DEPRECATED
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
                print "%.4f" % (math.exp(self.transitions[k][l])),
            print " "
         
    def getNumIP(self):
        """Return number of inheritance patterns"""
        return len(self.inheritance_patterns)
        
    def getNumPP(self):
        """Return number of phased inheritance patterns"""
        return len(self.phased_patterns)
        
    def getNumStates(self):
        """Return number of HMM states per one SNP position"""
        return len(self.states)   
    
    def isNormal(self, state):
        return state.inheritance_pattern == (1,1)
        
    def isReal(self, state):
        return isinstance(state.inheritance_pattern, tuple)
    
    def areRecombination(self, state1, state2):
        if state1.inheritance_pattern != state2.inheritance_pattern:
            return False
        eq = 0
        for hap in range(2):
            eq += (state1.phased_pattern[hap] == state2.phased_pattern[hap])
        return (eq > 0)
    
    def getCorrespondingInState(self, state):
        for st in self.states[self.getNumPP():]:
            if st.inheritance_pattern == state.inheritance_pattern and \
               st.phased_pattern == "in":
                return st
        return None
        
    def getCorrespondingOutState(self, state):
        for st in self.states[self.getNumPP():]:
            if st.inheritance_pattern == state.inheritance_pattern and \
               st.phased_pattern == "out":
                return st
        return None
        
    def getStartState(self):
        for i, state in enumerate(self.states):
            if state.inheritance_pattern == "s": return i, state
        return None
        
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
    
    def generateTransitionProb(self):
        """Generate transition probabilities between HMM states"""
        #number of possible states per one SNP position
        num_states = self.getNumStates()
        num_real_states = self.getNumPP() 
        
        #precompute statistics
        num_recombs = [0 for x in range(self.getNumIP())]
        num_cnvs = 0
        for state in self.states[:num_real_states]:
            num_recombs[self.inheritance_patterns.index(state.inheritance_pattern)] += 1
            num_cnvs += (not self.isNormal(state))
        
        
        trans = [[0. for x in range(num_states)] for y in range(num_states)]
        #first generate transitions from the real states
        for ip in self.inheritance_patterns:
            
            if ip == (1, 1): #generate Normal states transitions\
                pstay = 0.9499
                precomb = 0.05 / (num_recombs[self.inheritance_patterns.index(ip)] - 1)
                pgo = 0.0001
                for i, state1 in enumerate(self.states[:num_real_states]):
                    if state1.inheritance_pattern != ip: continue
                    #states "inside the IP component"
                    for j, state2 in enumerate(self.states[:num_real_states]):
                        if i == j: #stay in the state
                            trans[i][i] = pstay
                        elif state2.inheritance_pattern == ip: #to recombination of the same IP
                            trans[i][j] = precomb
                    #to the silent exit node    
                    outState = self.getCorrespondingOutState(state1)
                    trans[i][self.states.index(outState)] = pgo
                    
            else: #generate CNV states transitions
                pstay = 0.949
                precomb = 0.05 / (num_recombs[self.inheritance_patterns.index(ip)] - 1)
                pgo = 0.001
                for i, state1 in enumerate(self.states[:num_real_states]):
                    if state1.inheritance_pattern != ip: continue
                    #states "inside the IP component"
                    for j, state2 in enumerate(self.states[:num_real_states]):
                        if i == j: #stay in the state
                            trans[i][i] = pstay
                        elif state2.inheritance_pattern == ip: #to recombination of the same IP
                            trans[i][j] = precomb
                    #to the silent exit node    
                    outState = self.getCorrespondingOutState(state1)
                    trans[i][self.states.index(outState)] = pgo
        
        #now generate transitions from the silent states
        inNormal = self.getCorrespondingInState( HMMState((1,1), ()) )
        for i, state1 in enumerate(self.states[num_real_states:]):
            i += num_real_states
            #if it is the start node
            if state1.phased_pattern == "s":
                for j, state2 in enumerate(self.states[num_real_states:]):
                    if state2.phased_pattern == "in":
                        trans[i][num_real_states + j] = 1. / self.getNumIP()
            
            #if it is a silent exit node
            elif state1.phased_pattern == "out":
                exit_prob = 1. / self.getNumIP()
                j, st = self.getExitState()
                trans[i][j] = exit_prob
                        
                if self.isNormal(state1):
                    prob = 1. / self.getNumIP()
                    for j, state2 in enumerate(self.states[num_real_states:]):
                        if state2.phased_pattern == "in" and \
                         state2.inheritance_pattern != state1.inheritance_pattern:
                            trans[i][num_real_states + j] = prob
                else:
                    trans[i][self.states.index(inNormal)] = 1. - exit_prob
                    
            #if it is a silent starting node
            elif state1.phased_pattern == "in":
                prob = 1. / num_recombs[self.inheritance_patterns.index(state1.inheritance_pattern)]
                for j, state2 in enumerate(self.states[:num_real_states]):
                    if state2.inheritance_pattern == state1.inheritance_pattern:
                        trans[i][j] = prob
        
        #normalize and take logarithm
        for i in range(num_states):
            for j in range(num_states):
                if trans[i][j] < 10e-10: trans[i][j] = self.neg_inf
                else: trans[i][j] = math.log(trans[i][j])
            self.logNormalize(trans[i])
        
        return trans
    
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
        if val != None: return val
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
        if val != None: return val
        
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
        """
        Viterbi decoding of the most probable path.
        """
        num_states = self.getNumStates()
        num_real_states = self.getNumPP() 
        transitions = self.transitions
        
        n = len(samples)
        predecessor = [[0 for i in range(num_states)] for j in xrange(n+1)] 
        #DP table dim: seq length+1 x num_states
        table = [[self.neg_inf for i in range(num_states)] for j in xrange(n+1)] 
        
        start_id, start_state = self.getStartState()
        exit_id, exit_state = self.getExitState()
        
        
        '''INITIALIZE'''
        #the start state probability is 1 -> log(1)=0
        table[0][start_id] = 0. 
        #propagate over all silent states
        for state_id, state in enumerate(self.states[num_real_states:]):
            state_id += num_real_states
            if transitions[start_id][state_id] == self.neg_inf: continue #just for speed-up
            table[0][state_id] = table[0][start_id] + transitions[start_id][state_id]
            predecessor[0][state_id] = -1
         
        '''DYNAMIC PROGRAMMING'''
        #for all SNP positions do:
        for pos in xrange(1, n+1):
            #real positions are <1..n+1), pos 1 is the base case

            #(i) compute new values for all phased patterns - "real" states of the HMM:
            for state_id, state in enumerate(self.states[:num_real_states]):
                #emission probability in the given state
                emis_p = self.logLHGivenState(samples[pos-1], M[pos-1], P[pos-1], mixture, state)
                
                #porobability of this state is `emission * max{previous state * transition}`
                max_prev = self.neg_inf
                arg_max = -1
                for prev_id in range(num_states):
                    if transitions[prev_id][state_id] == self.neg_inf: continue #just for speed-up
                    tmp = table[pos-1][prev_id] + transitions[prev_id][state_id]
                    if tmp > max_prev:
                        max_prev = tmp
                        arg_max = prev_id
                table[pos][state_id] = max_prev + emis_p
                predecessor[pos][state_id] = arg_max
            
            #(ii, iii) transitions from 'real' states and silent states with lower id to silent states
            #note: the states in self.states are already ordered such that a transistion from a silent 
            #   state to another silent state always points to a state with heigher index
            for state_id, state in enumerate(self.states[num_real_states:]):
                state_id += num_real_states
                #porobability of this state is `max{(actual real state or silent state with lower id) * transition}`
                max_prev = self.neg_inf
                arg_max = -1
                for prev_id in range(state_id + 1):
                    if transitions[prev_id][state_id] == self.neg_inf: continue #just for speed-up
                    tmp = table[pos][prev_id] + transitions[prev_id][state_id]
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
        
        for i in xrange(1, n+1):
            state = self.states[ path[i] ]
            path[i] = self.inheritance_patterns.index(state.inheritance_pattern)
            
        return path[1:]
        
        
    def getP(self, v, u, j, i, samples, M, P, mixture):
        #base case
        if i == j: return float(v == u)
        
        #caching:
        code = (v, u, j)
        val = self.PCache.get(code)
        if val != None: return val
        
        num_real_states = self.getNumPP()
        transitions = self.transitions
        states = self.states
        logSum = self.logSum
        
        #for silent *v*
        if v >= num_real_states:
            summ = self.neg_inf
            for w in range(num_real_states):
                if transitions[v][w] == self.neg_inf: continue
                if states[v].inheritance_pattern != \
                   states[w].inheritance_pattern: continue
                tmp = transitions[v][w] + self.getP(w, u, j, i, samples, M, P, mixture)
                summ = logSum(summ, tmp)
            result = summ
            
        #for real *v*
        else:
            emis_p = self.logLHGivenState(samples[j-1], M[j-1], P[j-1], mixture, states[v])
            summ = self.neg_inf
            for w in range(num_real_states):
                if transitions[v][w] == self.neg_inf: continue
                if states[v].inheritance_pattern != \
                   states[w].inheritance_pattern: continue
                tmp = transitions[v][w] + self.getP(w, u, j+1, i, samples, M, P, mixture)
                summ = logSum(summ, tmp)
            result = summ + emis_p
        
        
        self.PCache[code] = result
        return result
    
    def preCompute(self, v, u, j, i, samples, M, P, mixture):
        code = (v, u, j)
        if i == j:
            self.PCache[code] = float(v == u)
        else:
            num_real_states = self.getNumPP()
            transitions = self.transitions
            states = self.states
            logSum = self.logSum
            
            #for silent *v*
            if v >= num_real_states:
                summ = self.neg_inf
                for w in range(num_real_states):
                    if transitions[v][w] == self.neg_inf: continue
                    if states[v].inheritance_pattern != \
                       states[w].inheritance_pattern: continue
                    tmp = transitions[v][w] + self.PCache.get((w, u, j))
                    summ = logSum(summ, tmp)
                result = summ
                
            #for real *v*
            else:
                emis_p = self.logLHGivenState(samples[j-1], M[j-1], P[j-1], mixture, states[v])
                summ = self.neg_inf
                for w in range(num_real_states):
                    if transitions[v][w] == self.neg_inf: continue
                    #if states[v].inheritance_pattern != \
                    #   states[w].inheritance_pattern: continue
                    tmp = transitions[v][w] + self.PCache.get((w, u, j+1))
                    summ = logSum(summ, tmp)
                result = summ + emis_p
            self.PCache[code] = result
    
    def extendedLabeling(self, samples, M, P, mixture):
        """
        Most probable extended labeling.
        """
        import sys
        print >> sys.stderr, "BEGIN"
        num_states = self.getNumStates()
        num_real_states = self.getNumPP() 
        transitions = self.transitions
        states = self.states
        getP = self.getP
        
        
        n = len(samples)
        predecessor = [[0 for i in range(num_states)] for j in xrange(n+1)]
        jump = [[0 for i in range(num_states)] for j in xrange(n+1)]  
        #DP table dim: seq length+1 x num_states
        table = [[self.neg_inf for i in range(num_states)] for j in xrange(n+1)] 
        
        start_id, start_state = self.getStartState()
        exit_id, exit_state = self.getExitState()
        
        
        '''INITIALIZE'''
        #the start state probability is 1 -> log(1)=0
        table[0][start_id] = 0. 
        #propagate over all silent states
        for state_id, state in enumerate(self.states[num_real_states:]):
            state_id += num_real_states
            if transitions[start_id][state_id] == self.neg_inf: continue #just for speed-up
            table[0][state_id] = table[0][start_id] + transitions[start_id][state_id]
            predecessor[0][state_id] = -1
            jump[0][state_id] = 0
         
        '''DYNAMIC PROGRAMMING'''
        #for all SNP positions do:
        for pos in xrange(1, n+1):
            #real positions are <1..n+1), pos 1 is the base case
            self.PCache = {}
            print >> sys.stderr, pos, 
            
            '''
            preCompute = self.preCompute
            for j in reversed(xrange(1,  pos+1)):
                i = pos
                for v in xrange(num_states):
                    for u in xrange(num_real_states):
                        preCompute(v, u, j, i, samples, M, P, mixture)
            '''            
            
            
            
            #(i) compute new values for all phased patterns - "real" states of the HMM:
            for state_id, state in enumerate(states[:num_real_states]):
                max_prev = self.neg_inf
                arg_max = -1
                jump_from = -1
                
                
                #if the pravious transition is the critical edge:
                #emission probability in current state
                emis_p = self.logLHGivenState(samples[pos-1], M[pos-1], P[pos-1], mixture, state)
                #porobability of this state is `emission * max{previous state * transition}`
                for prev_id in range(num_states):
                    if transitions[prev_id][state_id] == self.neg_inf: continue #just for speed-up
                    tmp = table[pos-1][prev_id] + transitions[prev_id][state_id] + emis_p
                    if tmp > max_prev:
                        max_prev = tmp
                        arg_max = prev_id
                        jump_from = pos-1             
                
                #where the previous critical edge is:
                for j in reversed(xrange(1,  pos)):
                    #between which states the critical edge is:
                    for to_id in range(num_real_states, num_states):
                        #the 'jumped' region after critical edge must have one lable
                        if states[to_id].inheritance_pattern != \
                           states[state_id].inheritance_pattern: continue
                        for from_id in range(num_real_states, num_states):
                            #there must be an edge
                            if transitions[from_id][to_id] == self.neg_inf: continue
                            #to be a critical edge, it must change the label
                            if states[from_id].inheritance_pattern == \
                               states[to_id].inheritance_pattern: continue
                            
                            #tmp = table[j-1][from_id] + transitions[from_id][to_id] + self.PCache.get((to_id, state_id, j))
                            tmp = table[j-1][from_id] + transitions[from_id][to_id] + getP(to_id, state_id, j, pos, samples, M, P, mixture)
                            if tmp > max_prev:
                                max_prev = tmp
                                arg_max = prev_id
                                jump_from = j-1
                
                #record the max value and traceback info
                table[pos][state_id] = max_prev
                predecessor[pos][state_id] = arg_max
                jump[pos][state_id] = jump_from
            
            #(ii, iii) transitions from 'real' states and silent states with lower id to silent states
            #note: the states in self.states are already ordered such that a transistion from a silent 
            #   state to another silent state always points to a state with heigher index
            for state_id, state in enumerate(states[num_real_states:]):
                state_id += num_real_states
                #porobability of this state is `max{(actual real state or silent state with lower id) * transition}`
                max_prev = self.neg_inf
                arg_max = -1
                jump_from = -1
                for prev_id in range(state_id + 1):
                    if transitions[prev_id][state_id] == self.neg_inf: continue #just for speed-up
                    tmp = table[pos][prev_id] + transitions[prev_id][state_id]
                    if tmp > max_prev:
                        max_prev = tmp
                        arg_max = prev_id
                        jump_from = pos
                table[pos][state_id] = max_prev
                predecessor[pos][state_id] = arg_max
                jump[pos][state_id] = jump_from
            
            #normalize to sum to 1
            #self.logNormalize(table[pos])
        
        '''RESULT'''
        print >> sys.stderr, "COMPUTING BACKTRACE"
        path = [-1 for i in xrange(n+1)]
        path[n] = exit_id 
        while path[n] >= num_real_states or path[n]!=start_id: path[n] = predecessor[n][path[n]]
        pos = n
        for x in range(jump[pos][path[pos]]+1, pos): path[x] = path[pos]
        while predecessor[pos][path[pos]] > -1:
            new_pos = jump[pos][path[pos]]
            path[new_pos] = predecessor[pos][path[pos]]
            pos = new_pos
            while path[pos] >= num_real_states or path[pos]!=start_id: path[pos] = predecessor[pos][path[pos]]
            for x in range(jump[pos][path[pos]]+1, pos): path[x] = path[pos]
        
        for i in xrange(1, n+1):
            state = states[ path[i] ]
            path[i] = self.inheritance_patterns.index(state.inheritance_pattern)
            
        return path[1:]
    
    def computeForward(self, samples, M, P, mixture):
        '''
        Posterior decoding: forward algorithm
        '''
        num_states = self.getNumStates()
        num_real_states = self.getNumPP()
        patterns = self.inheritance_patterns
        transitions = self.transitions

        n = len(samples)
        #scale factors
        scale = [0. for i in xrange(n+1)]
        #DP table dim: seq length+1 x num_states
        table = [[self.neg_inf for i in range(num_states)] for j in xrange(n+1)] 
        
        start_id, start_state = self.getStartState()
        exit_id, exit_state = self.getExitState()
        
        
        '''INITIALIZE'''
        #the start state probability is 1 -> log(1)=0
        table[0][start_id] = 0. 
        #propagate over all silent states
        for state_id, state in enumerate(self.states[num_real_states:]):
            state_id += num_real_states
            if transitions[start_id][state_id] == self.neg_inf: continue #just for speed-up
            table[0][state_id] = table[0][start_id] + transitions[start_id][state_id]
        
        '''DYNAMIC PROGRAMMING'''
        #for all SNP positions do:
        for pos in xrange(1, n+1):
            #real positions are <1..n+1), pos 1 is the base case
            
            #(i) compute new values for all phased patterns - "real" states of the HMM:
            for state_id, state in enumerate(self.states[:num_real_states]):
                #emission probability in the given state
                emis_p = self.logLHGivenState(samples[pos-1], M[pos-1], P[pos-1], mixture, state)
                
                #probability of this state is `emission * \sum {previous state * transition}`
                summ = self.neg_inf
                for prev_id in range(num_states):
                    if transitions[prev_id][state_id] == self.neg_inf: continue #just for speed-up
                    tmp =  table[pos-1][prev_id] + transitions[prev_id][state_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = emis_p + summ
            
            #(ii, iii) transitions from 'real' states and silent states with *lower* id to silent states
            #note: the states in self.states are already ordered such that a transistion from a silent 
            #   state to another silent state always points to a state with heigher index
            for state_id, state in enumerate(self.states[num_real_states:]):
                state_id += num_real_states
                
                #probability of this state is `\sum {(actual real state or silent state with lower id) * transition}`
                summ = self.neg_inf
                for prev_id in range(state_id):
                    if transitions[prev_id][state_id] == self.neg_inf: continue #just for speed-up
                    tmp =  table[pos][prev_id] + transitions[prev_id][state_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = summ
                
            #normalize to sum to 1
            scale[pos] = self.logNormalize(table[pos])
        
        '''RESULT'''
        #p(X) - probability of the observed sequence (samples)
        pX = table[n][exit_id]
        #print "fwd: ", math.exp(pX), pX

        return table, pX, scale
        
        
    def computeBackward(self, samples, M, P, mixture, scale):
        '''
        Posterior decoding: backward algorithm
        '''
        num_states = self.getNumStates()
        num_real_states = self.getNumPP()
        patterns = self.inheritance_patterns
        transitions = self.transitions

        n = len(samples)
        #DP table dim: seq length+2 x num_states
        table = [[self.neg_inf for i in range(num_states)] for j in xrange(n+1)]
        
        start_id, start_state = self.getStartState()
        exit_id, exit_state = self.getExitState()
        
        
        '''INITIALIZE'''
        #the exit state probability is 1 -> log(1)=0
        table[n][exit_id] = 0. 
        #(iii) add transitions from silent to silent states with *higher* id (*no emission*)
        for state_id in reversed(range(num_real_states, num_states)):
            #\sum {'next'_state * transition_from_here} (no emission)
            summ = self.neg_inf
            for next_id in reversed(range(state_id, num_states)):
                if transitions[state_id][next_id] == self.neg_inf: continue #just for speed-up
                tmp = transitions[state_id][next_id] + table[n][next_id]
                summ = self.logSum(summ, tmp)
            table[n][state_id] = self.logSum(summ, table[n][state_id])
        
        #(ii) add transitions from 'real' states to silent states (*no emission*)
        for state_id in range(num_real_states):
            #\sum {'next'_state * transition_from_here} (no emission)
            summ = self.neg_inf
            for next_id in range(num_real_states, num_states):
                if transitions[state_id][next_id] == self.neg_inf: continue #just for speed-up
                tmp = transitions[state_id][next_id] + table[n][next_id]
                summ = self.logSum(summ, tmp)
            table[n][state_id] = self.logSum(summ, table[n][state_id])
        
        '''DYNAMIC PROGRAMMING'''    
        #for all SNP positions do:
        for pos in reversed(xrange(0, n)):
            #real positions are <1..n+1), the 0 and n+1 are sentinels
            
            #precompute emission probabilities
            emis_p = []
            for state in self.states[:num_real_states]:
                emis_p.append(self.logLHGivenState(samples[pos], M[pos], P[pos], mixture, state) )
                
            #(i) transitions from all states to 'real' states of the HMM:
            for state_id in range(num_states):
                # \sum {emission_in_'next' * 'next'_state * transition_from_here}
                summ = self.neg_inf
                for next_id in range(num_real_states):
                    if transitions[state_id][next_id] == self.neg_inf: continue #just for speed-up
                    tmp = transitions[state_id][next_id] + table[pos+1][next_id] + emis_p[next_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = summ
            
            #(iii) add transitions from silent to silent states with *higher* id (*no emission*)
            for state_id in reversed(range(num_real_states, num_states)):
                #\sum {'next'_state * transition_from_here} (no emission)
                summ = self.neg_inf
                for next_id in reversed(range(state_id, num_states)):
                    if transitions[state_id][next_id] == self.neg_inf: continue #just for speed-up
                    tmp = transitions[state_id][next_id] + table[pos][next_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = self.logSum(summ, table[pos][state_id])
            
            #(ii) add transitions from 'real' states to silent states (*no emission*)
            for state_id in range(num_real_states):
                #\sum {'next'_state * transition_from_here} (no emission)
                summ = self.neg_inf
                for next_id in range(num_real_states, num_states):
                    if transitions[state_id][next_id] == self.neg_inf: continue #just for speed-up
                    tmp = transitions[state_id][next_id] + table[pos][next_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = self.logSum(summ, table[pos][state_id])
            
            #pseudonormalize by scaling factor from fwd
            for state_id in range(num_states):
                table[pos][state_id] -= scale[pos+1] 
        
        '''RESULT'''
        #p(X) - probability of the observed sequence (samples)
        pX = table[0][start_id]
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
            path.append( self.inheritance_patterns.index(max_IP) )
        
        return path
    
    
    def posteriorDecoding(self, samples, M, P, mixture):
        """Posterior probability decoding.
        
        For each position returns a list of states ordered by their posterior
        porobability. Returns list(list(tuple(probability, state)))
        
        """
        n = len(samples)
        num_real_states = self.getNumPP()
        
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
                max_ = self.neg_inf
                for s_id, s in enumerate(self.states[:num_real_states]):
                    if s.inheritance_pattern == ip:
                        #fwd, bck are in log space (therefore + instead *); p(X) is constant
                        max_ = max(max_, fwd[pos][s_id] + bck[pos][s_id])
                
                '''#sum over corresponding phased patterns
                sum_ = self.neg_inf
                for s_id, s in enumerate(self.states):
                    if s.inheritance_pattern == ip:
                        #fwd, bck are in log space (therefore + instead *); p(X) is constant
                        sum_ = self.logSum(sum_, fwd[pos][s_id] + bck[pos][s_id])
                '''
                    
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
        A = [[self.neg_inf for x in range(num_states)] for y in range(num_states)]
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
                    if self.transitions[k][l] == self.neg_inf:
                        A[k][l] = self.neg_inf
                        continue
                        
                    val = self.neg_inf
                    for pos in xrange(1, n):
                        emis_p = self.logLHGivenState(samples[pos], M[pos],\
                         P[pos], mixture, self.states[l])
                        tmp = fwd[pos][k] + self.transitions[k][l] + emis_p + bck[pos+1][l]
                        if tmp == self.neg_inf: continue
                        val = self.logSum(val, tmp)
                    
                    if val == self.neg_inf:
                        A[k][l] = self.neg_inf
                    else:
                        val -= pX1
                        A[k][l] = val
                    #if math.exp(A[k][l]) < 0.1:
                    #    A[k][l] = math.log(0.1)
                         
            #maximization
            old_trans = copy.deepcopy(self.transitions)
            for k in range(num_states):
                summ = self.neg_inf
                for m in range(num_states):
                    if A[k][m] == self.neg_inf: continue
                    summ = self.logSum(summ, A[k][m])
                #print summ, math.exp(summ)
                for l in range(num_states):
                    self.transitions[k][l] = A[k][l] - summ
                    #if self.transitions[k][l] < math.log(0.0005):
                    #    self.transitions[k][l] = math.log(0.0005)
            
            change = 0.
            for k in range(num_states):
                summ = self.neg_inf
                for m in range(num_states):
                    if A[k][m] == self.neg_inf: continue
                    summ = self.logSum(summ, self.transitions[k][m])
                #print summ
                for l in range(num_states):
                    #self.transitions[k][l] -= summ
                    if not (self.transitions[k][l] == self.neg_inf or old_trans[k][l] == self.neg_inf):
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
                        self.transitions[k][l] = self.neg_inf
                    else:
                        self.transitions[k][l] = math.log(A[k][l]) - summ
                    #if self.transitions[k][l] < math.log(0.0005):
                    #    self.transitions[k][l] = math.log(0.0005)
                        
            change = 0.
            for k in range(num_states):
                summ = self.neg_inf
                for m in range(num_states):
                    if A[k][m] == self.neg_inf: continue
                    summ = self.logSum(summ, self.transitions[k][m])
                for l in range(num_states):
                    #self.transitions[k][l] -= summ
                    if not (self.transitions[k][l] == self.neg_inf or old_trans[k][l] == self.neg_inf):
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
        num_real_states = self.getNumPP()
        
        table = []
        for pos in xrange(1, n+1):
            tmp = []
            for ip_id, ip in enumerate(self.inheritance_patterns):
                #max over corresponding phased patterns
                max_ = self.neg_inf
                for s in self.states[:num_real_states]:
                    if s.inheritance_pattern == ip:
                        ll = self.logLHGivenState(samples[pos-1], M[pos-1], P[pos-1], mixture, s)
                        max_ = max(max_, ll)
                
                '''#sum over corresponding phased patterns
                sum_ = self.neg_inf
                for s in self.states:
                    if s.inheritance_pattern == ip:
                        ll = self.logLHGivenState(samples[pos-1], M[pos-1], P[pos-1], mixture, s)
                        sum_ = self.logSum(sum_, ll)
                '''        

                tmp.append((max_ , ip_id))
                
            tmp.sort(reverse=True)    
            table.append(tmp)
        return table

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
