import random
import math
import numpypy as np
import itertools
import copy

nucleotides = ['A', 'C', 'G', 'T']

def readInput(file_samples_in, file_maternal_in, file_paternal_in):
    '''
    reads the input - maternal alleles, paternal alleles, and samples
    '''
    def parseHaplotypes(lines):
        X = []
        for line in lines:
            line = line.rstrip("\n")
            X.append(line.split(" "))
        return X
    def parseSamples(lines):
        X = []
        for line in lines:
            l = line.rstrip("\n").split(" ")
            for i in range(len(l)):
                l[i] = int(l[i])
            X.append(l)
        return X

    fpaternal = open(file_paternal_in, 'r')
    fmaternal = open(file_maternal_in, 'r')
    fsamples = open(file_samples_in, 'r')
    
    #maternal and paternal alleles
    M = parseHaplotypes(fmaternal.readlines())
    P = parseHaplotypes(fpaternal.readlines())
    
    #plasma samples
    samples = parseSamples(fsamples.readlines())
    
    return samples, M, P
    
class FCNV(object):
    def __init__(self):
        super(FCNV, self).__init__()
        
        #cache for log-likelihood values
        self.logLikelihoodCache = {}
        
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
        
        # trasition probabilities (uniform for all positions)
        p_stay = 0.95
        p_stay3 = 0.975
        self.transitions = [[math.log((1.-p_stay)/(self.num_states-0.)) for x in range(self.num_states)] for y in range(self.num_states)]
        for i in range(self.num_states):
            self.transitions[i][3] = math.log(2.*(1.-p_stay)/(self.num_states-0.))
        self.transitions[0][0] = math.log(p_stay)
        self.transitions[1][1] = math.log(p_stay)
        self.transitions[2][2] = math.log(p_stay)
        self.transitions[4][4] = math.log(p_stay)
        self.transitions[5][5] = math.log(p_stay)
        self.transitions[6][6] = math.log(p_stay)
        self.transitions[3] = [math.log((1.-p_stay3)/(self.num_states-1)) for x in range(self.num_states)]
        self.transitions[3][3] = math.log(p_stay3)
        
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
        dif = b-a
        #return a + math.log(1 + math.exp(b-a) )
        return a + precise(dif)
    
    def logMultinomial(self, xs, ps):
        '''
        [integers], [doubles] -> double
        calculates probability mass function of multinomial distribution
        returns log(Multi(xs | sum(xs), ps))
        >>> f = FCNV()
        >>> f.logMultinomial([1,2,3], [.2,.3,.5])
        -2.0024805005437063
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
            (NumPy array) -> (NumPy array)
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
        returns list of nucleotides probabilities at a locus with given
        maternal and fetal alleles, assuming the given mixture of fetal alleles
        
        >>> f = FCNV()
        >>> f.allelesDistribution(['A', 'A'], ['T', 'A'], 0.1)
        [0.92307692307692313, 0.0096153846153846159, 0.0096153846153846159, 0.057692307692307696]
        '''
        dist = []
        for nuc in nucleotides:
            temp = maternal.count(nuc) / float(len(maternal)) * (1.-mix)
            temp += fetal.count(nuc) / float(len(fetal)) * (mix)
            dist.append(temp)
        
        summ = 0.
        for x in dist:
            summ += 0.01 + x
        for i in range(len(dist)):
            dist[i] = (dist[i]+0.01)/summ #normalize
        return dist
            
    def logLikelihoodGivenPattern(self, nuc_counts, maternal_alleles, paternal_alleles, mix, pattern):
        '''
        >>> f = FCNV()
        >>> f.logLikelihoodGivenPattern([3, 0, 0, 71], ['T', 'T'], ['T', 'A'], 0.1, [1,1])
        -3.6974187244239021
        >>> f.logLikelihoodGivenPattern([3, 0, 0, 71], ['T', 'T'], ['T', 'A'], 0.1, [2,1])
        -3.4838423152150568
        >>> f.logLikelihoodGivenPattern([3, 0, 0, 71], ['T', 'T'], ['T', 'A'], 0.1, [0,2])
        -3.1614953504205019
        '''
        
        code = (tuple(nuc_counts), tuple(maternal_alleles), tuple(paternal_alleles), mix, tuple(pattern))
        val = self.logLikelihoodCache.get(code)
        if val: return val
        
        fetal_set = dict()
        set_size = 0
        for set_m in itertools.combinations(maternal_alleles, pattern[0]):
            for set_p in itertools.combinations(paternal_alleles, pattern[1]):
                if not fetal_set.get((set_m + set_p)):
                    fetal_set[(set_m + set_p)] = 1
                else:
                    fetal_set[(set_m + set_p)] +=1
                set_size += 1
                #print (set_m + set_p), fetal_set[(set_m + set_p)]
                
        result = "nothing"
        log_set_size = math.log(float(set_size))
        #print set_size
        for fetal_alleles in fetal_set:
                ps = self.allelesDistribution(maternal_alleles, fetal_alleles, mix)
                tmp = math.log(fetal_set[fetal_alleles]) - log_set_size + self.logMultinomial(nuc_counts, ps)
                result = self.logSum(result, tmp)
                #print "%0.5f" %multinomial(nuc_counts, ps), nuc_counts, \
                #maternal_alleles, fetal_alleles,  ["%0.4f" % i for i in ps]
        
        self.logLikelihoodCache[code] = result
        return result
    
    def estimateMixture(self, samples, M, P):
        '''
        estimate mixture of fetal genome in the plasma samples
        '''
        mix = 0.
        num = 0
        for i in xrange(len(samples)):
            if M[i][0] == M[i][1] and not P[i][0] == P[i][1]:
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
            emis_p = self.logLikelihoodGivenPattern(samples[0], M[0], P[0], mixture, patterns[i])
            table[1][i] = math.log(1./num_states)+ emis_p #initial state probabilities   
            predecessor[(1, i)] = -1
        
        #for all SNP positions do:
        for pos in xrange(2, n+1):
            #real positions are <1..n+1), pos 1 is the base case

            #for each SNP consider all inheritance patterns - states of the HMM:
            for state in range(num_states):
                #emission probability in the given state
                emis_p = self.logLikelihoodGivenPattern(samples[pos-1], M[pos-1], P[pos-1],\
                 mixture, patterns[state])
                #print emis_p
                
                #porobability of this state is emission * max{previous state * transition}
                max_prev = float("-inf")
                arg_max = -1
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
        table = [[0. for i in range(num_states)] for j in xrange(n+1)] 

        #initital base case
        for i in range(num_states):
            emis_p = self.logLikelihoodGivenPattern(samples[0], M[0], P[0],\
                 mixture, patterns[i])
            table[1][i] = math.log(1./num_states)+ emis_p #initial state probabilities

        #for all SNP positions do:
        for pos in xrange(2, n+1):
            #real positions are <1..n+1), pos 1 is the base case

            #for each SNP consider all inheritance patterns - states of the HMM:
            for state in range(num_states):
                #emission probability in the given state
                emis_p = self.logLikelihoodGivenPattern(samples[pos-1], M[pos-1], P[pos-1],\
                 mixture, patterns[state])
                #print emis_p
                
                #probability of this state is emission * \sum {previous state * transition}
                summ = "nothing"
                for s in range(num_states):
                    tmp = transitions[s][state]+table[pos-1][s]
                    summ = self.logSum(summ, tmp)
                table[pos][state] = emis_p + summ
        
        #p(X) - probability of the observed sequence (samples)
        pX = "nothing"
        for i in table[n]:
            pX = self.logSum(pX, i)
        print "fwd: ", np.exp(pX), pX

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
        table = [[0. for i in range(num_states)] for j in xrange(n+1)]
        
        #initital base case
        for i in range(num_states):
            table[n][i] = 0. #initial state probabilities are 1 -> log(1)=0
            
        #for all SNP positions do:
        for pos in reversed(xrange(n)):
            #real positions are <1..n+1), the 0 and n+1 are sentinels
            if pos == 0: continue
            
            #precompute emission probabilities
            emis_p = []
            for s in range(num_states):
                emis_p.append(self.logLikelihoodGivenPattern(samples[pos], M[pos], P[pos], mixture, patterns[s]) )
            
            #for each SNP consider all inheritance patterns - states of the HMM:
            for state in range(num_states):
                #probability of this state is  \sum {emission *previous state * transition_to_here}
                summ = "nothing"
                for s in range(num_states):
                    tmp = transitions[state][s] + table[pos+1][s]+ emis_p[s]
                    summ = self.logSum(summ, tmp)
                table[pos][state] = summ
        
        #p(X) - probability of the observed sequence (samples)
        pX = "nothing"
        for s in range(num_states):
            emis_p = self.logLikelihoodGivenPattern(samples[0], M[0], P[0], mixture, patterns[s])
            tmp = math.log(1./num_states) + table[1][s] + emis_p
            pX = self.logSum(pX, tmp)
        #print "bck: ", np.exp(pX), pX
        
        return table, pX
        
    
    def maxPosteriorDecoding(self, samples, M, P, mixture):
        '''
        Maximum posterior probability decoding
        '''
        n = len(samples)
        
        fwd, pX1 = self.computeForward(samples, M, P, mixture)
        bck, pX2 = self.computeBackward(samples, M, P, mixture)
        
        if abs(pX1-pX2) > 0.1:
            print "p(X) from forward and backward DP doesn't match", pX1, pX2
            return
        
        path = []
        for pos in xrange(1, n+1):
            max_posterior = float("-inf")
            arg_max = -1
            for s in range(self.num_states):
                # fwd, bck ate in log space (therefore + instead *); p(X) is constant
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
        
        if abs(pX1-pX2) > 0.1:
            print "p(X) from forward and backward DP doesn't match", pX1, pX2
            return
        
        path = []
        for pos in xrange(1, n+1):
            tmp = []
            for s in range(self.num_states):
                # fwd, bck ate in log space (therefore + instead *); p(X) is constant
                tmp.append((fwd[pos][s] + bck[pos][s], s))
            path.append(reversed(sorted(tmp)))
        return path
        
    def baumWelshForTransitions(self, samples, M, P, mixture):
        '''
        Baum-Welsh EM algorithm for HMM parameter estimation,
        but implemented only for traning transition probabilities
        '''
        n = len(samples)
        A = [[0. for x in range(self.num_states)] for y in range(self.num_states)]
        num_iter = 0
        while True and num_iter < 3:
            num_iter += 1
            #expectation
            fwd, pX1 = self.computeForward(samples, M, P, mixture)
            bck, pX2 = self.computeBackward(samples, M, P, mixture)
            if abs(pX1-pX2) > 0.1:
                print "p(X) from forward and backward DP doesn't match", pX1, pX2
                return  
            
            
            for k in range(self.num_states):
                for l in range(self.num_states):
                    val = "nothing"
                    for pos in xrange(1, n):
                        emis_p = self.logLikelihoodGivenPattern(samples[pos], M[pos],\
                         P[pos], mixture, self.inheritance_patterns[l])
                        tmp = fwd[pos][k] + self.transitions[k][l] + emis_p + bck[pos+1][l]
                        val = self.logSum(val, tmp)
                    val -= pX1
                    A[k][l] = val
                    #if np.exp(A[k][l]) < 0.1:
                    #    A[k][l] = math.log(0.1)
                         
            #maximization
            change = 0.
            for k in range(self.num_states):
                summ = "nothing"
                for m in range(self.num_states):
                    summ = self.logSum(summ, A[k][m])
                #print summ, np.exp(summ)
                for l in range(self.num_states):
                    change += (np.exp(max(self.transitions[k][l], math.log(0.001)))\
                     - np.exp(A[k][l] - summ))**2
                    self.transitions[k][l] = A[k][l] - summ
                    if np.exp(self.transitions[k][l]) < 0.001:
                        self.transitions[k][l] = math.log(0.001)
            print "change: ", change 
            if change < 0.01: break
                    
        tt = np.array(self.transitions)
        print np.ceil(np.exp(tt)*1000)/1000.
        #for k in xrange(self.num_states):
        #    for l in xrange(self.num_states):
        #        print np.exp(self.transitions[k][l])
        #    print " "
            
    def viterbiTrainingForTransitions(self, samples, M, P, mixture):
        '''
        Viterbi training algorithm for HMM parameter estimation,
        but implemented only for traning transition probabilities
        '''
        n = len(samples)
        A = [[0. for x in range(self.num_states)] for y in range(self.num_states)]
        for opak in range(2):
            #expectation
            path, table = self.viterbiPath(samples, M, P, mixture)
            fwd, pX1 = self.computeForward(samples, M, P, mixture)
            print pX1
            
            for k in range(self.num_states):
                for l in range(self.num_states):
                    val = 0
                    for pos in xrange(1, n-1):
                        if path[pos] == k and path[pos+1] == l:
                            val += 1
                    if val < 3: val = 3
                    A[k][l] = val
                    
                         
            #maximization
            for k in range(self.num_states):
                summ = math.log(sum(A[k]))
                #for m in range(self.num_states):
                #    summ = self.logSum(summ, A[k][m])
                #print summ, np.exp(summ)
                for l in range(self.num_states):
                    self.transitions[k][l] = math.log(A[k][l]) - summ
                    if np.exp(self.transitions[k][l]) < 0.001:
                        self.transitions[k][l] = math.log(0.001)
                    
        tt = np.array(self.transitions)
        print np.ceil(np.exp(tt)*1000)/1000.
        #for k in xrange(self.num_states):
        #    for l in xrange(self.num_states):
        #        print np.exp(self.transitions[k][l])
        #    print " "     

def test(fcnv, samples, M, P, mixture, file_name):
    vp, vt = fcnv.viterbiPath(samples, M, P, mixture)
    posterior = fcnv.posteriorDecoding(samples, M, P, mixture)
    
    ground_truth = []
    fetal_in = open(file_name, 'r')
    for line in fetal_in.readlines():
        line = line.strip("\n")
        ground_truth.append(line.split(" ")[-1])
        
    print fcnv.inheritance_patterns
    viterbi_correct = 0
    max_posterior_correct = 0
    posterior_dist = [0,0,0,0,0,0,0]
    state_stats = [0,0,0,0,0,0,0]
    state_correct_vp = [0,0,0,0,0,0,0]
    state_correct_pp = [0,0,0,0,0,0,0]
    for i in xrange(len(vp)):
        state_stats[int(ground_truth[i])] += 1
        post = []
        for x in posterior[i]:   
            post.append(x[1])
            
        #print ground_truth[i], vp[i], pp[i], '|', post
        viterbi_correct += int(int(ground_truth[i]) == vp[i])
        state_correct_vp[int(ground_truth[i])] += int(int(ground_truth[i]) == vp[i])
        max_posterior_correct += int(int(ground_truth[i]) == post[0])
        state_correct_pp[int(ground_truth[i])] += int(int(ground_truth[i]) == post[0])
        
        x = post.index(int(ground_truth[i]))
        posterior_dist[x] += 1
        #print ground_truth[i], vp[i], post, '|', post.index(int(ground_truth[i]))
        
    print 'Viterbi  : ',(viterbi_correct*100.)/len(vp), '% OK'
    print 'Posterior: ',(max_posterior_correct*100.)/len(vp), '% OK'
    posterior_dist = np.array(posterior_dist)
    posterior_dist = (posterior_dist*100.)/len(vp)
    print posterior_dist
    
    print "stats  : ", state_stats
    print "viterbi: ", state_correct_vp
    print "mposter: ", state_correct_pp

def main():
    mixture = 0.1 #proportion of fetal genome in plasma
    file_fetal_out = "fetal_prediction.txt"    
    #ffetal = open(file_fetal_out, 'w')

    #read the input
    samples_t, M_t, P_t = readInput("plasma_samplesT.txt", "maternal_allelesT.txt", "paternal_allelesT.txt")
    samples, M, P = readInput("plasma_samples.txt", "maternal_alleles.txt", "paternal_alleles.txt")
    fcnv = FCNV()
    
    mix_t = fcnv.estimateMixture(samples_t, M_t, P_t)
    mix = fcnv.estimateMixture(samples, M, P)
    
    #mix = mix_t = 0.1 #proportion of fetal genome in plasma
    print "Est. Mixture: ", mix_t, mix
    #return
    
    print "---------- BEFORE TRAINING --------------"
    print ">>>> ON TRAINING DATA:"
    test(fcnv, samples_t, M_t, P_t, mix_t, "fetal_allelesT.txt")
    print ">>>> ON TESTING DATA:"
    test(fcnv, samples, M, P, mix, "fetal_alleles.txt")
    print "========================================="

    def_trans = copy.deepcopy(fcnv.transitions)
    print "---------- AFTER up to 3 rounds of B-W --------------"
    fcnv.baumWelshForTransitions(samples_t, M_t, P_t, mix_t)
    print ">>>> ON TRAINING DATA:"
    test(fcnv, samples_t, M_t, P_t, mix_t, "fetal_allelesT.txt")
    print ">>>> ON TESTING DATA:"
    test(fcnv, samples, M, P, mix, "fetal_alleles.txt")
    print "========================================="
    
    fcnv.transitions = copy.deepcopy(def_trans)
    print "------- AFTER 2 rounds of Viterbi Training ----------"
    fcnv.viterbiTrainingForTransitions(samples_t, M_t, P_t, mix_t)
    print ">>>> ON TRAINING DATA:"
    test(fcnv, samples_t, M_t, P_t, mix_t, "fetal_allelesT.txt")
    print ">>>> ON TESTING DATA:"
    test(fcnv, samples, M, P, mix, "fetal_alleles.txt")
    print "========================================="
    
    
    #fcnv.baumWelshForTransitions(samples, M, P, mixture)
    #fcnv.viterbiTrainingForTransitions(samples, M, P, mixture)
    #print fcnv.computeForward(samples, M, P, mixture)
    #print fcnv.computeBackward(samples, M, P, mixture)
    #pp = fcnv.maxPosteriorDecoding(samples, M, P, mixture)

    #test(fcnv, samples, M, P, mixture, "fetal_alleles.txt")
    
    '''
    n=74+1
    summ = "nothing"
    for x1 in xrange(0, n):
      for x2 in xrange(0, n-x1):
        for x3 in xrange(0, n-x1-x2):
          for x4 in xrange(0, n-x1-x2-x3):
            if x1+x2+x3+x4 == n-1:
              tmp = fcnv.logLikelihoodGivenPattern([x1, x2, x3, x4], ['T', 'T'], ['T', 'A'], 0.15, [1,2])
              summ = fcnv.logSum(summ, tmp)
    print summ, np.exp(summ)
    '''
    
if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    
    main()
    



