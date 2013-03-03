import math
import numpypy as np
import sys
import fcnvHMM

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

def computeEval(reference, prediction, sensitivity):
    num_real = 0
    num_called = 0
    ok_called = []
    not_called = []
    
    ind = 0
    while ind < len(reference):
        ctype = reference[ind]
        length = 0
        begin_pos = ind
        while ind < len(reference) and reference[ind] == ctype:
            length += 1
            ind += 1
        end_pos = ind
        
        max_ok_length = 0
        tmp_length = 0
        stats = [0, 0, 0, 0, 0, 0, 0] 
        for i in range(begin_pos, end_pos):
            stats[prediction[i]] += 1
            if prediction[i] == ctype:
                tmp_length += 1
            else:
                #print i, prediction[i],  ctype, type(prediction[i]), type(ctype)
                max_ok_length = max(max_ok_length, tmp_length)
                tmp_length = 0
        max_ok_length = max(max_ok_length, tmp_length)
        
        num_real += 1
        if max_ok_length >= length/float(sensitivity):
            num_called += 1
            ok_called.append((length, ctype))
        else:
            max_predic = -1
            max_type = -1
            for i in range(len(stats)):
                if stats[i] > max_predic:
                    max_predic = stats[i]
                    max_type = i
            not_called.append((length, ctype, max_type))
        #else: print  max_ok_length, length, reference[ind], prediction[ind], ctype
    return num_called, num_real, ok_called, not_called

def categorize(data, by_length, by_type):
    lengths = [100, 500, 1000, 5000, 10000, 20000, 50000, 100000]
    num_patterns = 7
    for item in data:
        length = item[0]
        ptype = item[1]
        gen_len = sorted([(abs(length-l), l) for l in lengths])[0][1]
        #print length, gen_len
        if not gen_len in by_length: by_length[gen_len] = 0
        by_length[gen_len] += 1
        if not ptype in by_type: by_type[ptype] = 0
        by_type[ptype] += 1

def computeStats(ok, wrong, pref):
    lengths = [100, 500, 1000, 5000, 10000, 20000, 50000, 100000]
    ok_by_length = dict([(i,0) for i in lengths])
    ok_by_type = dict([(i,0) for i in range(7)])
    wr_by_length = dict([(i,0) for i in lengths])
    wr_by_type = dict([(i,0) for i in range(7)])
    categorize(ok, ok_by_length, ok_by_type)
    categorize(wrong, wr_by_length, wr_by_type)
    for l in lengths:
        o = ok_by_length[l]
        w = wr_by_length[l]
        ratio = 0
        if o+w > 0: ratio = o/float(o+w) *100
        print pref, l, ": ", o, w, ratio, '%'
    for t in range(7):
        o = ok_by_type[t]
        w = wr_by_type[t]
        ratio = 0
        if o+w > 0: ratio = o/float(o+w) *100
        print pref, t, ": ", o, w, ratio, '%'

def test(fcnv, samples, M, P, mixture, ground_truth):
    vp, vt = fcnv.viterbiPath(samples, M, P, mixture)
    posterior = fcnv.posteriorDecoding(samples, M, P, mixture)
    posterior = [list(x) for x in posterior]

    byLL = fcnv.likelihoodDecoding(samples, M, P, mixture)    

    print fcnv.inheritance_patterns
    
    viterbi_correct = 0
    max_posterior_correct = 0
    ll_correct = 0
    #states distributions for posterior and pure likelihood
    posterior_dist = [0,0,0,0,0,0,0]
    ll_dist = [0,0,0,0,0,0,0]
    #other stats
    state_stats = [0,0,0,0,0,0,0]
    state_correct_vp = [0,0,0,0,0,0,0]
    state_correct_pp = [0,0,0,0,0,0,0]
    state_correct_ll = [0,0,0,0,0,0,0]
    pp = []
    tableLH = []
    avg_distance = 0.
    distance_set = {}
    ll_misclass = [[0]*7 for x in range(7)]
    state3_stats = [{ (0,0):0, (0,1):0, (1,0):0, (1,1):0 } for x in range(7)]
    for i in xrange(len(vp)):
        state_stats[ground_truth[i]] += 1
        post = []
        for x in posterior[i]:   
            post.append(x[1])
        pp.append(post[0])
        
        ll_state = []
        ll_value = []
        for x in byLL[i]:   
            ll_state.append(x[1])
            ll_value.append(x[0])
        tableLH.append(ll_value)
            
        #print ground_truth[i], vp[i], pp[i], '|', post
        viterbi_correct += int(ground_truth[i] == vp[i])
        state_correct_vp[ground_truth[i]] += int(ground_truth[i] == vp[i])
        max_posterior_correct += int(ground_truth[i] == post[0])
        state_correct_pp[ground_truth[i]] += int(ground_truth[i] == post[0])
        x = post.index(ground_truth[i])
        posterior_dist[x] += 1
        
        #for pure likelihood
        
        threshold = .1
        x = ll_state.index(ground_truth[i])
        ll_diff = abs(ll_value[x]-ll_value[0])
        if ll_diff < threshold:
            ll_dist[0] += 1
            #ll_correct += int(ground_truth[i] == ll_state[0])
            ll_correct += 1
            #state_correct_ll[ground_truth[i]] += int(ground_truth[i] == ll_state[0])
            state_correct_ll[ground_truth[i]] += 1
        else:            
            ll_dist[x] += 1
            for j in range(7):
                if ll_value[x] - ll_value[j] < 0.01:
                #if abs(ll_value[0] - ll_value[j]) < 1.:
                    ll_misclass[ground_truth[i]][ll_state[j]] += 1
                    if ground_truth[i] == 3:
                        m_poly = M[i][0] == M[i][1]
                        p_poly = P[i][0] == P[i][1]
                        state3_stats[ll_state[j]][(m_poly, p_poly)] += 1
        
        avg_distance += ll_diff
        int_distance = int(round(ll_diff))
        if int_distance not in distance_set:
            distance_set[int_distance] = 1
        else:
            distance_set[int_distance] += 1
        #print ground_truth[i], vp[i], post, '|', post.index(int(ground_truth[i]))
    
    '''
    for i in xrange(len(vp)):
        print ground_truth[i], vp[i],  
        emis = []
        for j in range(fcnv.num_states):
            emis_p = fcnv.logLikelihoodGivenPattern(samples[i], M[i], P[i], mixture, fcnv.inheritance_patterns[j])
            emis.append((emis_p, j))
        print samples[i],
        for e in reversed(sorted(emis)):
            print str(e[1])+":%0.1f" % e[0] ,
        print " "
        print "   ",
        for x in posterior[i]:
            print str(x[1])+":%0.1f" % x[0] ,
        print " "
        #print tableLH[i]
    '''    
        
    #np.set_printoptions(precision=3)
    posterior_dist = np.array(posterior_dist)
    posterior_dist = (posterior_dist*100.)/len(vp)    
    print posterior_dist  
    ll_dist = np.array(ll_dist)
    ll_dist = (ll_dist*100.)/len(vp)    
    print ll_dist
    print "stats of misclassified by LL:"
    ll_misclass = np.array(ll_misclass, float)
    for i in range(len(ll_misclass)):
        ll_misclass[i] = 100.*ll_misclass[i]/float(state_stats[i])
        for j in range(len(ll_misclass[i])): ll_misclass[i][j] = round(ll_misclass[i][j]*100)/100.
    print ll_misclass
    print "mistakes for state3: "
    for x in range(7):
        print x,":", state3_stats[x]
    print "stats  : ", state_stats
    print "viterbi: ", state_correct_vp
    print "mposter: ", state_correct_pp
    print "byLHood: ", state_correct_ll
    print "avgDist: ", avg_distance / len(vp) 
    print "dist_set:", distance_set
    
    
    print 'Viterbi  : ',(viterbi_correct*100.)/len(vp), '% OK'
    print 'Posterior: ',(max_posterior_correct*100.)/len(vp), '% OK'
    print 'LikeliH. : ',(ll_correct*100.)/len(vp), '% OK'
    
    for i in [2]: #, 2, 4, 8, 16]:
        #precision and recall
        print "sensitivity: 1 /", i
        
        #recall Viterbi
        called_v, real, ok_called, not_called = computeEval(ground_truth, vp, i)
        print "viterbi recall   : ", called_v, '/',  real, " => ", 
        print "%0.3f" % (called_v*100./real), "%"
        print ok_called
        print not_called
        computeStats(ok_called, not_called, "VR")
        
        #recall max posterior
        called_p, real, ok_called, not_called = computeEval(ground_truth, pp, i)       
        print "mposter recall   : ", called_p, '/',  real, " => ",
        print "%0.3f" % (called_p*100./real), "%"
        print ok_called
        print not_called
        computeStats(ok_called, not_called, "PR")
        
        #prediction Viterbi 
        correct_v, claimed_v, ok_prediction, wr_prediction = computeEval(vp, ground_truth, i)
        print "viterbi precision: ", correct_v , '/',  claimed_v, " => ", 
        print "%0.3f" % (correct_v*100./claimed_v), "%"
        print ok_prediction
        print wr_prediction
        computeStats(ok_prediction, wr_prediction, "VP")
        
        #prediction max posterior
        correct_p, claimed_p, ok_prediction, wr_prediction = computeEval(pp, ground_truth, i)
        print "mposter precision: ", correct_p , '/',  claimed_p, " => ", 
        print "%0.3f" % (correct_p*100./claimed_p), "%"
        print ok_prediction
        print wr_prediction
        computeStats(ok_prediction, wr_prediction, "PP")
        
        #format for LaTeX tabular
        #print "%0.3f" % ((viterbi_correct*100.)/len(vp)), "\%"
        #print "%0.3f" % ((max_posterior_correct*100.)/len(vp)), "\%"
        print "%0.3f" % (called_v*100./real), "\%"
        print "%0.3f" % (called_p*100./real), "\%"
        print "%0.3f" % (correct_v*100./claimed_v), "\%"
        print "%0.3f" % (correct_p*100./claimed_p), "\%"
        
        
def main():
    mixture = 0.1 #proportion of fetal genome in plasma
    #file_fetal_out = "fetal_prediction.txt"    
    #ffetal = open(file_fetal_out, 'w')
    
    #parse command line arguments => suffix of the output files
    path = sys.argv[0]
    path = path[:path.rfind('/')+1]
    path_in = path + "in/"
    suffix = ""
    if len(sys.argv) == 2: suffix = sys.argv[1]
    
    #read the input
    '''
    samples_t, M_t, P_t = readInput( \
        path_in + "plasma_samplesT" + suffix + ".txt", \
        path_in + "maternal_allelesT" + suffix + ".txt", \
        path_in + "paternal_allelesT" + suffix + ".txt")
    '''
    samples, M, P = readInput( \
        path_in + "plasma_samples" + suffix + ".txt", \
        path_in + "maternal_alleles" + suffix + ".txt", \
        path_in + "paternal_alleles" + suffix + ".txt")

    fcnv = fcnvHMM.FCNV()
    
    #mix_t = fcnv.estimateMixture(samples_t, M_t, P_t)
    mix = fcnv.estimateMixture(samples, M, P)
    
    #mix = mix_t = 0.15 #proportion of fetal genome in plasma
    #print "Est. Mixture: ", mix_t, mix
    print "Est. Mixture: ", mix
    
    '''
    ground_truth_t = []
    fetal_in = open(path_in + "fetal_allelesT" + suffix + ".txt", 'r')
    for line in fetal_in.readlines():
        line = line.rstrip("\n")
        ground_truth_t.append(int(line.split(" ")[-1]))
    fetal_in.close()
    '''
    
    ground_truth = []
    fetal_in = open(path_in + "fetal_alleles" + suffix + ".txt", 'r')
    for line in fetal_in.readlines():
        line = line.rstrip("\n")
        ground_truth.append(int(line.split(" ")[-1]))
    fetal_in.close()
    
    print "------------------ w/o TRAINING -------------------"
    test(fcnv, samples, M, P, mix, ground_truth)
    
    '''
    print "------------------ BEFORE TRAINING -------------------"
    print ">>>>>>>>>>>>>>>>>>>>>>>>> ON TRAINING DATA: >>>>>>>>>>>>>>>>>>>>>"
    test(fcnv, samples_t, M_t, P_t, mix_t, ground_truth_t)
    print ">>>>>>>>>>>>>>>>>>>>>>>>> ON TESTING DATA: >>>>>>>>>>>>>>>>>>>>>>"
    test(fcnv, samples, M, P, mix, ground_truth)
    print "=================================================================================="
    
    
    def_trans = copy.deepcopy(fcnv.transitions)
    print "----------------- AFTER up to 5 rounds of B-W -------------------"
    fcnv.baumWelshForTransitions(samples_t, M_t, P_t, mix_t)
    print ">>>>>>>>>>>>>>>>>>>>>>>>> ON TRAINING DATA: >>>>>>>>>>>>>>>>>>>>>"
    test(fcnv, samples_t, M_t, P_t, mix_t, ground_truth_t)
    print ">>>>>>>>>>>>>>>>>>>>>>>>> ON TESTING DATA: >>>>>>>>>>>>>>>>>>>>>>"
    test(fcnv, samples, M, P, mix, ground_truth)
    print "=================================================================================="
    
    fcnv.transitions = copy.deepcopy(def_trans)
    print "----------- AFTER up to 5 rounds of Viterbi Training ------------"
    fcnv.viterbiTrainingForTransitions(samples_t, M_t, P_t, mix_t)
    print ">>>>>>>>>>>>>>>>>>>>>>>>> ON TRAINING DATA: >>>>>>>>>>>>>>>>>>>>>"
    test(fcnv, samples_t, M_t, P_t, mix_t, ground_truth_t)
    print ">>>>>>>>>>>>>>>>>>>>>>>>> ON TESTING DATA: >>>>>>>>>>>>>>>>>>>>>>"
    test(fcnv, samples, M, P, mix, ground_truth)
    print "=================================================================================="
    '''
    
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
    print summ, math.exp(summ)
    '''
    
if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    
    main()
    



