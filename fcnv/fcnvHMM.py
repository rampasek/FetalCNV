import random
import math
import numpypy as np
import itertools
import copy

nucleotides = ['A', 'C', 'G', 'T']

class FCNV(object):
    def __init__(self):
        super(FCNV, self).__init__()
        
        #cache for log-likelihood values
        self.logLikelihoodCache = {}
        self.distributionCache = {}
        
        #precompute lookup table for logSum function        
        self.logTable1 = []
        for i in range(501):
            self.logTable1.append(math.log(1 + math.exp(-i/100.)) )
        self.logTable2 = []
        for i in range(51):
            self.logTable2.append(math.log(1 + math.exp(-i)) )
        
        #generate inheritance patterns
        self.inheritance_patterns = []
        for mats in range(3):
            for pats in range(3):
                if mats+pats in [0,4]: continue
                self.inheritance_patterns.append([mats, pats])
        
        # number of possible states per one SNP position
        self.num_states = len(self.inheritance_patterns) 
        
        # fdsalfjdasl
        self.history_depth = 3
        
        # trasition probabilities (uniform for all positions)
        p_stay = 0.999
        p_stay3 = 0.9999
        #self.transitions = [[math.log((1.-p_stay)/(self.num_states-0.)) for x in range(self.num_states)] for y in range(self.num_states)]
        self.transitions = [[float("-inf") for x in range(self.num_states)] for y in range(self.num_states)]
        for i in range(self.num_states):
            self.transitions[i][3] = math.log(1.-p_stay)
        self.transitions[0][0] = math.log(p_stay)
        self.transitions[1][1] = math.log(p_stay)
        self.transitions[2][2] = math.log(p_stay)
        self.transitions[4][4] = math.log(p_stay)
        self.transitions[5][5] = math.log(p_stay)
        self.transitions[6][6] = math.log(p_stay)
        self.transitions[3] = [math.log((1.-p_stay3)/(self.num_states-1)) for x in range(self.num_states)]
        self.transitions[3][3] = math.log(p_stay3)
        
        self.main_state = 3
        #print the table    
        #for k in xrange(self.num_states):
        #    for l in xrange(self.num_states):
        #        print "%0.4f" % (math.exp(self.transitions[k][l])), 
        #    print " "
            
        
    def logSum(self, x, y):
        '''
        returns sum of two numbers in log space
        i.e. equivalent to log(exp(x)+exp(y))
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
        
        if x == "nothing": return y
        elif y == "nothing": return x
       
        a = max(x, y)
        b = min(x, y)
        #dif = b-a
        return a + math.log(1 + math.exp(b-a) )
    
    def logMultinomial(self, xs, ps):
        '''
        [integers], [doubles] -> double
        calculates probability mass function of multinomial distribution
        returns log(Multi(xs | sum(xs), ps))
        >>> f = FCNV()
        >>> f.logMultinomial([1,2,3], [.2,.3,.5])
        -2.0024805005437072
        '''
        
        def gammaln(n):
            '''
            Logarithm of Euler's gamma function for discrete values.
            '''
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
            '''
            (float) -> (float)
            returns ln(x!)
            '''
            if isinstance(x, list):
                res = []
                for val in x:
                    res.append(gammaln(val+1))
                return res
            else: 
                return gammaln(x+1)
        
        n = sum(xs)
        #xs, ps = np.array(xs), np.array(ps)
        #result = logFactorial(n) - sum(logFactorial(xs)) + sum(xs * np.log(ps))
        
        result = logFactorial(n) - sum(logFactorial(xs))
        for i in range(len(ps)):
            result += xs[i] * math.log(ps[i])
        return result

    def allelesDistribution(self, maternal, fetal, mix):
        '''
        returns list of nucleotides probabilities at a position with given
        maternal and fetal alleles, assuming the given mixture of fetal alleles
        
        >>> f = FCNV()
        >>> f.allelesDistribution(['A', 'A'], ['T', 'A'], 0.1)
        [0.9230769230769231, 0.009615384615384616, 0.009615384615384616, 0.057692307692307696]
        '''
        #code = (tuple(maternal), tuple(fetal), mix)
        #print code
        #val = self.distributionCache.get(code)
        #if val: return val
        
        adjusted_fetal_admix = mix/2. * len(fetal)
        adjusted_maternal_admix = (1.-mix)/2. * len(maternal)
        cmix = adjusted_fetal_admix / (adjusted_maternal_admix + adjusted_fetal_admix)
        dist = []
        for nuc in nucleotides:
            p = maternal.count(nuc) / float(len(maternal)) * (1.-cmix)
            p += fetal.count(nuc) / float(len(fetal)) * (cmix)
            dist.append(p)
        
        summ = 0.
        for x in dist:
            summ += 0.01 + x
        for i in range(len(dist)):
            dist[i] = (dist[i]+0.01)/summ #normalize
        
        #self.distributionCache[code] = dist
        return dist
            
    def logLikelihoodGivenPhasedPattern(self, nuc_counts, maternal_alleles, fetal_alleles, mix):
        '''
        >>> f = FCNV()
        >>> f.logLikelihoodGivenPattern([3, 0, 0, 71], ['T', 'T'], ['T', 'A'], 0.1, [1,1])
        -3.6974187244239025
        >>> f.logLikelihoodGivenPattern([3, 0, 0, 71], ['T', 'T'], ['T', 'A'], 0.1, [2,1]) #-3.4838423152150568
        -3.6507523341244412
        >>> f.logLikelihoodGivenPattern([3, 0, 0, 71], ['T', 'T'], ['T', 'A'], 0.1, [0,2])
        -3.1614953504205028
        '''
        
        code = (tuple(nuc_counts), tuple(maternal_alleles), tuple(fetal_alleles), mix)
        #print code
        val = self.logLikelihoodCache.get(code)
        if val: return val
        
        ps = self.allelesDistribution(maternal_alleles, fetal_alleles, mix)
        result = self.logMultinomial(nuc_counts, ps)
        #print "%0.5f" %self.logMultinomial(nuc_counts, ps), nuc_counts, \
        #maternal_alleles, fetal_alleles,  ["%0.4f" % i for i in ps], fetal_set[fetal_alleles]
        
        self.logLikelihoodCache[code] = result
        return result
        
    def logLikelihoodGivenPattern(self, samples, M, P, pos, mix, pattern):
        '''
        smooth out likelihood by taking max average over previous positions
        '''
        #go over all maternal and paternal haplotype patterns possible for given inheritance pattern
        history_depth = self.history_depth
        max_ll = float("-inf")
        likelihood_sum = "nothing"
        num_configs = 0;
        for mHpatt in itertools.combinations_with_replacement([0,1], pattern[0]):
            for pHpatt in itertools.combinations_with_replacement([0,1], pattern[1]):
                likelihood = 0. #"nothing"
                num_configs += 1
                for i in range(pos-history_depth+1, pos+1):
                    if i < 0: continue
                    if i >= len(M): continue
                    fetal_alleles = []
                    for x in mHpatt: fetal_alleles.append(M[i][x])
                    for x in pHpatt: fetal_alleles.append(P[i][x])
                    fetal_alleles = tuple(sorted(fetal_alleles))
                    ll = self.logLikelihoodGivenPhasedPattern(samples[i], M[i], fetal_alleles, mix)
                    likelihood += ll #self.logSum(likelihood, ll)
                max_ll = max(max_ll, likelihood) # - math.log(history_depth))
                likelihood_sum = self.logSum(likelihood_sum, likelihood) # - math.log(history_depth))
        return likelihood_sum - math.log(num_configs)
        '''
        likelihood_sum = "nothing"
        max_ll = float("-inf")
        configs = 0
        for mHpatt in itertools.combinations_with_replacement([0,1], pattern[0]):
            for pHpatt in itertools.combinations_with_replacement([0,1], pattern[1]):
                fetal_alleles = []
                for x in mHpatt: fetal_alleles.append(maternals[x])
                for x in pHpatt: fetal_alleles.append(paternals[x])
                fetal_alleles = tuple(sorted(fetal_alleles))
                ll = self.logLikelihoodGivenPhasedPattern(counts, maternals, fetal_alleles, mix)
                
                likelihood_sum = self.logSum(likelihood_sum, ll)
                configs += 1
                max_ll = max(max_ll, ll)
        return max_ll #likelihood_sum - math.log(configs)
        '''
        
    def estimateMixture(self, samples, M, P):
        '''
        estimate mixture of fetal genome in the plasma samples
        '''
        mix = 0.
        num = 0
        for i in xrange(len(samples)):
            if M[i][0] == M[i][1] and not P[i][0] == P[i][1] and \
              (P[i][0] == M[i][0] or P[i][1] == M[i][0]):
                #print M[i], P[i], samples[i]            
                type1 = M[i][0]
                type2 = P[i][0]
                if type1 == type2: type2 = P[i][1]
                num_type2 = samples[i][nucleotides.index(type2)]

                if not num_type2 == 0:
                    sum_all = sum(samples[i])
                    num += 1
                    mix += (2.*num_type2)/sum_all
        return mix/num
    
    def viterbiPath(self, samples, M, P, mixture):
        '''
        Viterbi decoding of the most probable path
        '''
        num_states = self.num_states
        patterns = self.inheritance_patterns
        transitions = self.transitions
        
        n = len(samples)
        predecessor = {}
        #DP table dim: seq length+1 x num_states
        table = [[0. for i in range(num_states)] for j in xrange(n+1)] 
        
        #initital base case
        for i in range(num_states):
            emis_p = self.logLikelihoodGivenPattern(samples, M, P, 0, mixture, patterns[i])
            table[1][i] = math.log(1./num_states)+ emis_p #initial state probabilities   
            predecessor[(1, i)] = -1
        
        #for all SNP positions do:
        for pos in xrange(2, n+1):
            #real positions are <1..n+1), pos 1 is the base case

            #for each SNP consider all inheritance patterns - states of the HMM:
            for state in range(num_states):
                #emission probability in the given state
                emis_p = self.logLikelihoodGivenPattern(samples, M, P, pos-1, mixture, patterns[state])
                #print emis_p
                
                #porobability of this state is emission * max{previous state * transition}
                max_prev = float("-inf")
                arg_max = -1
                if not state == self.main_state:
                    max_prev = table[pos-1][self.main_state] + transitions[self.main_state][state]
                    arg_max = self.main_state
                    tmp = table[pos-1][state] + transitions[state][state]
                    if tmp >= max_prev:
                        max_prev = tmp
                        arg_max = state
                else:
                    for s in range(num_states):
                        tmp = table[pos-1][s] + transitions[s][state]
                        if tmp > max_prev:
                            max_prev = tmp
                            arg_max = s
                table[pos][state] = max_prev + emis_p
                predecessor[(pos, state)] = arg_max
        
        path = [-1 for i in xrange(n+1)]
        maxx = max(table[n])
        path[n] = table[n].index(maxx)
        for i in reversed(xrange(n)):
            path[i] = predecessor[(i+1, path[i+1])]
        
        '''
        print patterns
        for i in range(len(path)-1):
            i = i+1
            print "%2d" %(i), ': ', patterns[path[i]], path[i]
        '''   
        return path[1:], table
    
    def computeForward(self, samples, M, P, mixture):
        '''
        Posterior decoding: forward algorithm
        '''
        num_states = self.num_states
        patterns = self.inheritance_patterns
        transitions = self.transitions
        
        n = len(samples)
        #DP table dim: seq length+1 x num_states
        table = [[float("-inf") for i in range(num_states)] for j in xrange(n+1)] 

        #initital base case
        for i in range(num_states):
            emis_p = self.logLikelihoodGivenPattern(samples, M, P, 0, mixture, patterns[i])
            table[1][i] = math.log(1./num_states)+ emis_p #initial state probabilities

        #for all SNP positions do:
        for pos in xrange(2, n+1):
            #real positions are <1..n+1), pos 1 is the base case

            #for each SNP consider all inheritance patterns - states of the HMM:
            for state in range(num_states):
                #emission probability in the given state
                emis_p = self.logLikelihoodGivenPattern(samples, M, P, pos-1, mixture, patterns[state])
                #print emis_p
                
                #probability of this state is emission * \sum {previous state * transition}
                summ = "nothing"
                if not state == self.main_state:
                    summ = table[pos-1][self.main_state] + transitions[self.main_state][state]
                    tmp = table[pos-1][state] + transitions[state][state]
                    summ = self.logSum(summ, tmp)
                else:
                    for s in range(num_states):
                        tmp =  table[pos-1][s] + transitions[s][state]
                        summ = self.logSum(summ, tmp)
                table[pos][state] = emis_p + summ
        
        #p(X) - probability of the observed sequence (samples)
        pX = "nothing"
        for i in table[n]:
            pX = self.logSum(pX, i)
        print "fwd: ", math.exp(pX), pX

        return table, pX
        
    def computeBackward(self, samples, M, P, mixture):
        '''
        Posterior decoding: backward algorithm
        '''
        num_states = self.num_states
        patterns = self.inheritance_patterns
        transitions = self.transitions
        
        n = len(samples)
        #DP table dim: seq length+2 x num_states
        table = [[float("-inf") for i in range(num_states)] for j in xrange(n+1)]
        
        #initital base case
        for i in range(num_states):
            table[n][i] = 0. #initial state probabilities are 1 -> log(1)=0
            
        #for all SNP positions do:
        for pos in reversed(range(1, n)):
            #real positions are <1..n+1), the 0 and n+1 are sentinels
            
            #precompute emission probabilities
            emis_p = []
            for s in range(num_states):
                emis_p.append(self.logLikelihoodGivenPattern(samples, M, P, pos, mixture, patterns[s]) )
            
            #for each SNP consider all inheritance patterns - states of the HMM:
            for state in range(num_states):
                #probability of this state is  \sum {emission *previous state * transition_to_here}
                summ = "nothing"
                if not state == self.main_state:
                    summ = transitions[state][self.main_state] + table[pos+1][self.main_state] + emis_p[self.main_state]
                    tmp = transitions[state][state] + table[pos+1][state] + emis_p[state]
                    summ = self.logSum(summ, tmp)
                else:
                    for s in range(num_states):
                        tmp = transitions[state][s] + table[pos+1][s] + emis_p[s]
                        summ = self.logSum(summ, tmp)
                table[pos][state] = summ
        
        #p(X) - probability of the observed sequence (samples)
        pX = "nothing"
        for s in range(num_states):
            emis_p = self.logLikelihoodGivenPattern(samples, M, P, 0, mixture, patterns[s])
            tmp = math.log(1./num_states) + table[1][s] + emis_p
            pX = self.logSum(pX, tmp)
        #print "bck: ", math.exp(pX), pX
        
        return table, pX
        
    
    def maxPosteriorDecoding(self, samples, M, P, mixture):
        '''
        Maximum posterior probability decoding
        '''
        n = len(samples)
        
        fwd, pX1 = self.computeForward(samples, M, P, mixture)
        bck, pX2 = self.computeBackward(samples, M, P, mixture)
        
        if abs(pX1-pX2) > 0.01:
            print "p(X) from forward and backward DP doesn't match", pX1, pX2
            return
        
        path = []
        for pos in xrange(1, n+1):
            max_posterior = float("-inf")
            arg_max = -1
            for s in range(self.num_states):
                # fwd, bck are in log space (therefore + instead *); p(X) is constant
                tmp = fwd[pos][s] + bck[pos][s]
                if tmp > max_posterior:
                    max_posterior = tmp
                    arg_max = s
            #print max_posterior
            path.append(arg_max)
        
        return path
  
    def posteriorDecoding(self, samples, M, P, mixture):
        '''
        posterior probability decoding, for each position returns a list of states
        ordered by their posterior porobability
        '''
        n = len(samples)
        
        fwd, pX1 = self.computeForward(samples, M, P, mixture)
        bck, pX2 = self.computeBackward(samples, M, P, mixture)
        
        if abs(pX1-pX2) > 0.01:
            print "p(X) from forward and backward DP doesn't match", pX1, pX2
            return
        
        table = []
        for pos in xrange(1, n+1):
            tmp = []
            for s in range(self.num_states):
                # fwd, bck are in log space (therefore + instead *); p(X) is constant
                tmp.append((fwd[pos][s] + bck[pos][s], s))
            table.append(reversed(sorted(tmp)))
        return table
        
    def baumWelshForTransitions(self, samples, M, P, mixture):
        '''
        Baum-Welsh EM algorithm for HMM parameter estimation,
        but implemented only for traning transition probabilities
        '''
        n = len(samples)
        A = [[float("-inf") for x in range(self.num_states)] for y in range(self.num_states)]
        num_iter = 0
        while True and num_iter < 5:
            num_iter += 1
            #expectation
            fwd, pX1 = self.computeForward(samples, M, P, mixture)
            bck, pX2 = self.computeBackward(samples, M, P, mixture)
            if abs(pX1-pX2) > 0.01:
                print "p(X) from forward and backward DP doesn't match", pX1, pX2
                return  
            
            
            for k in range(self.num_states):            
                for l in range(self.num_states):
                    if self.transitions[k][l] == float("-inf"):
                        A[k][l] = float("-inf")
                        continue
                        
                    val = "nothing"
                    for pos in xrange(1, n):
                        emis_p = self.logLikelihoodGivenPattern(samples, M, P, pos, \
                            mixture, self.inheritance_patterns[l])
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
            for k in range(self.num_states):
                summ = "nothing"
                for m in range(self.num_states):
                    if A[k][m] == float("-inf"): continue
                    summ = self.logSum(summ, A[k][m])
                #print summ, math.exp(summ)
                for l in range(self.num_states):
                    self.transitions[k][l] = A[k][l] - summ
                    #if self.transitions[k][l] < math.log(0.0005):
                    #    self.transitions[k][l] = math.log(0.0005)
            
            change = 0.
            for k in range(self.num_states):
                summ = "nothing"
                for m in range(self.num_states):
                    if A[k][m] == float("-inf"): continue
                    summ = self.logSum(summ, self.transitions[k][m])
                #print summ
                for l in range(self.num_states):
                    #self.transitions[k][l] -= summ
                    if not (self.transitions[k][l] == float("-inf") or old_trans[k][l] == float("-inf")):
                        change += (self.transitions[k][l] - old_trans[k][l])**2
                    
            print "change: ", change 
            if change < 1: break
                    
        #print the table    
        for k in xrange(self.num_states):
            for l in xrange(self.num_states):
                print "%0.4f" % (math.exp(self.transitions[k][l])), 
            print " "
            
    def viterbiTrainingForTransitions(self, samples, M, P, mixture):
        '''
        Viterbi training algorithm for HMM parameter estimation,
        but implemented only for traning transition probabilities
        '''
        n = len(samples)
        A = [[0. for x in range(self.num_states)] for y in range(self.num_states)]
        num_iter = 0
        while True and num_iter < 5:
            num_iter += 1
            #expectation
            path, table = self.viterbiPath(samples, M, P, mixture)
            p_path = max(table[len(samples)])
            print "P of the max probable path: ", p_path
            
            for k in range(self.num_states):
                for l in range(self.num_states):
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
            for k in range(self.num_states):
                summ = sum(A[k])
                if summ == 0:
                    summ = float("inf")
                else: 
                    summ = math.log(sum(A[k]))
                #for m in range(self.num_states):
                #    summ = self.logSum(summ, A[k][m])
                #print summ, math.exp(summ)
                for l in range(self.num_states):
                    if A[k][l] == 0:
                        self.transitions[k][l] = float("-inf")
                    else:
                        self.transitions[k][l] = math.log(A[k][l]) - summ
                    #if self.transitions[k][l] < math.log(0.0005):
                    #    self.transitions[k][l] = math.log(0.0005)
                        
            change = 0.
            for k in range(self.num_states):
                summ = "nothing"
                for m in range(self.num_states):
                    if A[k][m] == float("-inf"): continue
                    summ = self.logSum(summ, self.transitions[k][m])
                for l in range(self.num_states):
                    #self.transitions[k][l] -= summ
                    if not (self.transitions[k][l] == float("-inf") or old_trans[k][l] == float("-inf")):
                        change += (self.transitions[k][l] - old_trans[k][l])**2
                    
            print "change: ", change 
            if change < 1: break
                   
        #print the table    
        for k in xrange(self.num_states):
            for l in xrange(self.num_states):
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
            for s in range(self.num_states):
                ll = self.logLikelihoodGivenPattern(samples, M, P, pos-1, \
                    mixture, self.inheritance_patterns[s])
                tmp.append((ll, s))
            table.append(reversed(sorted(tmp)))
        return table


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
