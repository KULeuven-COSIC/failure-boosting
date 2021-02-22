import numpy as np
from dist import Dist
from scipy.stats import norm
from math import log, pi, cos, sin, sqrt, acos
import tqdm
import itertools
from geometric import Ptheta, PcosthetaSurvival, PthetaMax, listToAlphabeta
from scipy.interpolate import RegularGridInterpolator

ACCURACY = 2**-300
# grid 1st failure
GRIDPOINTS1 = 100
GRIDPOINTSS1 = 30

# grid
GRIDPOINTS = 100
GRIDPOINTTHETACE = 100

# precalculation time
GRIDPOINTSPRECALC = 100
GRIDPOINTSPRECALC2 = 100
GRIDPOINTST1 = 100
GRIDPOINTSTHETA = 100

GRIDPOINTSS = 50
GRIDPOINTSEPP = 20

def calcFailgivennormctNoInfo_uneven(normc1, normc2, dist_norms1, dist_norms2, eprimeprime, thres, n2, pcosthetaSurvival, dist_costheta):
    """ Calculate the failure probability given a certain norm of c1 and c2.
    Takes into account distribution norm of S and e'' distribution.
    Approximations: geometric approximation
    :param normc: norm of C
    :param normsdist: distribution of norm of S
    :param epp: distribution of e''
    :param thres: failure threshold (typically q/4)
    :param n2: coefficients in transmitted messages
    :param pcosthetaSurvival: survival function of theta
    """
    # loop over all values and average failure probabilities
    p_fail = 0
    for norms1 in dist_norms1.keys():
        for norms2 in dist_norms2.keys():
            sc2 = normc2 * norms2
            p_tmp = 0
            for costheta1, p_costheta1 in dist_costheta.items():
                sc1 = norms1 * normc1 * costheta1
                for g, p_g in eprimeprime.items():
                    thres1 = (thres - g - sc1) / sc2
                    thres2 = (thres + g + sc1) / sc2
                    p_tmp += p_g * p_costheta1 * (pcosthetaSurvival(thres1) +  pcosthetaSurvival(thres2))
            p_fail += p_tmp * dist_norms1[norms1] * dist_norms2[norms2]

    # take into account n2 coefficients
    if p_fail > 2**-10:
        p_fail = 1 - (1 - p_fail)**n2
    else:
        p_fail = p_fail * n2
    return p_fail


def firstgeometric_uneven(n, n2, thres, s, sprime, e, eprime, eprimeprime):
    """ calculate alpha-beta list for first failure using geometric approximation. More accurate method to use is gaussian approximation. This implementation is for reference purposes.
    Approximations: 
    - geometric approximation, 
    - equal distribution of s and e and equal distribution of s' and e', 
    - approximates some distribution to GRIDPOINTS values for efficiency reasons, 
    - does not failure boost with e'' (it is taken into account when calculating failure probability but not for determining if a ciphertext is weak) 
    :param n: dimension
    :param n2: coefficients in transmitted message
    :param thres: failure threshold (typically q/4)
    :param s: distribution of s
    :param sprime: distribution of s'
    :param e: distribution of e
    :param eprime: distribution of e'
    :param eprimeprime: distribution of e''
    """
    ########################################################################
    # Make s, e, s', e' symmetric so geometric approximation holds
    # Add asymetric part to e'' and a possible constant part to thres (this will result in a small underestimation of the efficiency of failure boosting)
    ########################################################################
    s, mean_s = s.symmetry_split()
    sprime, mean_sprime = sprime.symmetry_split()
    e, mean_e = e.symmetry_split()
    eprime, mean_eprime = eprime.symmetry_split()
    
    eprimeprime += (mean_s * eprime + mean_eprime * s - mean_e * sprime - mean_sprime * e)**n

    thres += n * mean_s.mean() * mean_eprime.mean() 
    thres -= n * mean_e.mean() * mean_sprime.mean() 

    ########################################################################
    # Make a grid of possible normc1, normc2 combinations and calculate the probability
    ########################################################################

    # make a grid of possible ciphertext points
    dist_normc1 = sprime.norm(n).simplify(GRIDPOINTS1)
    xrange1 = dist_normc1.keys()

    dist_normc2 = eprime.norm(n).simplify(GRIDPOINTS1)
    xrange2 = dist_normc2.keys()

    # get the probability of ciphertexts over the grid
    p_cipher = []
    for normc1, normc2 in tqdm.tqdm(itertools.product(xrange1, xrange2), total=len(xrange1) * len(xrange2), desc='P[ciphertext]: '):
        prob = dist_normc1[normc1] * dist_normc2[normc2]
        p_cipher.append(prob)

    # normalize the probability
    sumprob = sum(p_cipher)
    p_cipher = [i / sumprob for i in p_cipher]

    ########################################################################
    # calculate the failure probability of all ciphertexts in the normc1, normc2 list
    ########################################################################

    # calculate norms
    dist_norms1 = e.norm(n).simplify(GRIDPOINTSS1).normalize()
    dist_norms2 = s.norm(n).simplify(GRIDPOINTSS1).normalize()

    # get angle distribution
    dist_costheta = Dist({ float(cos(theta)) : Ptheta(theta, n) for theta in np.linspace(0, float(pi), 100000) })
    dist_costheta = dist_costheta.normalize().simplify(GRIDPOINTSTHETA).normalize()

    # simplify epp distribution
    eprimeprime = eprimeprime.simplify(GRIDPOINTSEPP).normalize()

    # make a list of ciphertext probability and failure probability
    pcosthetaSurvival = PcosthetaSurvival(n)
    thelist = []
    for i, tmp in tqdm.tqdm(enumerate(itertools.product(xrange1, xrange2)), total=len(xrange1) * len(xrange2), desc='make list: '):
        normc1, normc2 = tmp
        if p_cipher[i] > ACCURACY:
            fail = calcFailgivennormctNoInfo_uneven(normc1, normc2, dist_norms1, dist_norms2, eprimeprime, thres, n2, pcosthetaSurvival, dist_costheta)
            thelist.append((p_cipher[i], fail))

    ########################################################################
    # convert the list to alpha beta values
    ########################################################################
    alpha, beta = listToAlphabeta(thelist)

    return alpha, beta


def normSgivenft_uneven(dist_norms1_apriori, dist_norms2_apriori, mean_normc1_0, mean_normc2_0, epp, n, n2, thres, targets, pcosthetaSurvival):
    dist_costheta = Dist({ float(cos(theta)) : Ptheta(theta, n) for theta in np.linspace(0, float(pi), 100000) })
    dist_costheta = dist_costheta.normalize().simplify(GRIDPOINTSTHETA)

    # determine P[F | ft] and P[F | norms, beta0]
    PFgivenSft = Dist({})
    for norms1 in tqdm.tqdm(dist_norms1_apriori.keys()):
        PFgivenSft[norms1] = calcFailgivennormctNoInfo_uneven(mean_normc1_0, mean_normc2_0, {norms1 : 1.}, dist_norms2_apriori, epp, thres, n2, pcosthetaSurvival, dist_costheta)

    PFgivenSft = PFgivenSft.normalize()
    PFgivenft = sum([PFgivenSft[norms]*dist_norms1_apriori[norms] for norms in dist_norms1_apriori])

    # calculate new normsdist
    normsdistNew = {}
    for norms1 in dist_norms1_apriori.keys():
        normsdistNew[norms1] = dist_norms1_apriori[norms1] * PFgivenSft[norms1] / (PFgivenSft[norms1] + (targets-1) * PFgivenft)

    normsdistNew = Dist(normsdistNew).normalize()

    return normsdistNew


def moregeometric_uneven(n, n2, thres, s, sprime, e, eprime, eprimeprime, n_failure, targets=1, beta0=0):
    """calculate directional failure boosting using geometric approach"""

    ########################################################################
    # Make s, e, s', e' symmetric so geometric approximation holds
    # Add asymetric part to e'' and a possible constant part to thres (this will result in a small underestimation of the efficiency of failure boosting)
    ########################################################################
    s, mean_s = s.symmetry_split()
    sprime, mean_sprime = sprime.symmetry_split()
    e, mean_e = e.symmetry_split()
    eprime, mean_eprime = eprime.symmetry_split()
    
    eprimeprime += (mean_s * eprime + mean_eprime * s - mean_e * sprime - mean_sprime * e)**n

    thres += n * mean_s.mean() * mean_eprime.mean() 
    thres -= n * mean_e.mean() * mean_sprime.mean() 

    eprimeprime = eprimeprime.simplify(GRIDPOINTSEPP)

    # distribution might not behave well when taking power or norm, whithout simplification calculation might take very long
    e.set_simplify_level(2**12)
    s.set_simplify_level(2**12)
    sprime.set_simplify_level(2**12)
    eprime.set_simplify_level(2**12)


    ########################################################################
    # Make a grid of possible normc1, normc2, theta_CE1, theta_CE2 combinations
    ########################################################################
    dist_normc1 = sprime.norm(n).simplify(GRIDPOINTS)
    xrange1 = dist_normc1.keys()
    
    dist_normc2 = eprime.norm(n).simplify(GRIDPOINTS)
    xrange2 = dist_normc2.keys()

    dist_theta = Dist({ angle : Ptheta(angle, n) for angle in np.linspace(0, float(pi), 10000) })
    dist_theta = dist_theta.normalize()
    dist_thetaSE = dist_theta.simplify(GRIDPOINTTHETACE)
    yrange = dist_thetaSE.keys()

    ########################################################################
    # Determine the mean norm of c1, c0 in first stage
    ########################################################################

    # calculate normc1_0 and normc2_0 value due to failure boosting 1
    if beta0 == None:
        mean_normc1_0 = dist_normc1.mean()
        mean_normc2_0 = dist_normc2.mean()
    else:
        e.set_simplify_level(2**12)
        s.set_simplify_level(2**12)
        means1 = e.norm(n).mean()
        means2 = s.norm(n).mean()
        list_tmp = []
        # calculate a list of normc1_0 normc2_0 values with increasing failure probability
        for normc1_0, normc2_0 in itertools.product(dist_normc1, dist_normc2):
            z = normc1_0 * means1 + normc2_0 * means2
            p_cipher = dist_normc1[normc1_0] * dist_normc2[normc2_0]
            list_tmp.append((z, p_cipher, normc1_0, normc2_0))

        # select the subset (of probability beta) of normc1_0, normc2_0 with highest failure probabity
        list_tmp.sort(key=lambda x: x[0], reverse=True)
        beta = 0
        normclist = []
        for z, p_cipher, normc1_0, normc2_0 in list_tmp:
            normclist.append((p_cipher, normc1_0, normc2_0))
            beta += p_cipher
            if beta > beta0:
                break

        # calculate the mean of the normc1_0, normc2_0 given beta
        normalizer = sum( [p_cipher for p_cipher, normc1, normc2 in normclist] )
        mean_normc1_0 = sum( [p_cipher * normc1 for p_cipher, normc1, normc2 in normclist] ) / normalizer
        mean_normc2_0 = sum( [p_cipher * normc2 for p_cipher, normc1, normc2 in normclist] ) / normalizer

    ########################################################################
    # Take into account weak secret if multiple targets are attacked in first stage
    ########################################################################

    if targets>1:
        pcosthetaSurvival = PcosthetaSurvival(n)
        dist_norms1_apriori = e.norm(n).simplify(GRIDPOINTS)
        dist_norms2_apriori = s.norm(n).simplify(GRIDPOINTS)
        dist_norms1 = normSgivenft_uneven(dist_norms1_apriori, dist_norms2_apriori, mean_normc1_0, mean_normc2_0, eprimeprime, n, n2, thres, targets, pcosthetaSurvival)
        dist_norms2 = normSgivenft_uneven(dist_norms2_apriori, dist_norms1_apriori, mean_normc2_0, mean_normc1_0, eprimeprime, n, n2, thres, targets, pcosthetaSurvival)
    else:
        dist_norms1 = e.norm(n).simplify(GRIDPOINTSS)
        dist_norms2 = s.norm(n).simplify(GRIDPOINTSS)
    
    ########################################################################
    # Determine theta_SE
    ########################################################################
    meanNorms1 = dist_norms1.mean()
    meanNorms2 = dist_norms2.mean()

    thetaCS = acos( thres / (mean_normc1_0 * meanNorms1 + mean_normc2_0 * meanNorms2) )
    costhetaSE = cos(thetaCS) / sqrt( cos(thetaCS)**2 + sin(thetaCS)**2 / (n_failure-1) )
    thetaSE = float( acos(costhetaSE) )
    pcosthetaSurvival = PcosthetaSurvival(2*n - 1)

    ########################################################################
    # calculate the probability of all ciphertexts in the grid
    ########################################################################
    # for more accurate results one can take into account that there are n2 coefficients and thus there is a higher probability of finding high values of cos(theta1) and cos(theta2) when considering multiple positions and thus angles. Not done due to efficiency reasons. The impact should be limited (but feel free to try for yourself)
    p_cipher = []
    cosdist_thetaSE = {}
    for theta, prob in dist_theta.simplify(1000).items():
        cosdist_thetaSE[cos(theta)] = cosdist_thetaSE.get(cos(theta), 0.) + prob
    cosdist_thetaSE = Dist(cosdist_thetaSE)

    for normc1, normc2 in tqdm.tqdm(itertools.product(xrange1, xrange2), total=len(xrange1) * len(xrange2), desc='P[ciphertext]: '):
        # zterm1 = (normc1 * meanNorms1 * cos(thetaCS))
        # zterm2 = (normc2 * meanNorms2 * cos(thetaCS))
        # zdist = cosdist_thetaSE * zterm1 + cosdist_thetaSE * zterm2
        # zmaxdist = zdist.simplify(10000).cf()

        sumprob_z = 0
        p_cipher_tmp = []
        for theta1, theta2 in itertools.product(yrange, yrange):
            # z = cos(theta1) * zterm1 + cos(theta2) * zterm2
            p_z = dist_thetaSE[theta1] * dist_thetaSE[theta2] #* zmaxdist(z)**(n2-1)
            prob =  dist_normc1[normc1] * dist_normc2[normc2] * p_z
            sumprob_z += p_z 
            p_cipher_tmp.append(prob)
        p_cipher_tmp = [ i / sumprob_z for i in p_cipher_tmp ]
        p_cipher.extend(p_cipher_tmp)

    sumprob = sum(p_cipher)

    # normalize the probability
    p_cipher = [i / sumprob for i in p_cipher]

    ########################################################################
    # calculate error probability if |cos(theta_j)| < |cos(theta_i)| 
    # this is the residual error probability of ciphertexts at different locations
    ########################################################################

    # This step is only done in the geometric approach, not here due to efficiency reasons. However, the impact should be limited as especially when directional failure boosting the error probability of |cos(theta_i)| will greatly and increasingly be higher than the one for the other angles. Maximal error is log2(n2) bits, but in practice it will be negligible

    ########################################################################
    # precalculate the failure probability of a1 * cos(theta1) + a2 * cos(theta2) > b + epp
    # for given a1, a2, b and random theta1, theta2, epp
    ########################################################################

    # make grid
    a1_max = dist_normc1.max() * dist_norms1.max() * sin(thetaCS)
    a2_max = dist_normc2.max() * dist_norms2.max() * sin(thetaCS)

    a1_values = np.linspace(0, a1_max+1, GRIDPOINTSPRECALC)
    a2_values = np.linspace(0.1, a2_max+1, GRIDPOINTSPRECALC)
    b_values = np.linspace( - 2 * thres, 2 * thres, GRIDPOINTSPRECALC2)
    V = np.zeros((GRIDPOINTSPRECALC, GRIDPOINTSPRECALC, GRIDPOINTSPRECALC2))

    # preparation functions
    cost1_dist = Dist({cos(t) : Ptheta(t, n) for t in np.linspace(0, pi, 1000)})
    cost1_dist = cost1_dist.normalize().simplify(GRIDPOINTST1)
    pcosthetaSurvival = PcosthetaSurvival(n)

    # loop over all a1, a2, b values and determine failure probability
    for idx1, a1 in tqdm.tqdm(enumerate(a1_values), total = GRIDPOINTSPRECALC, desc= 'precalc failure'):
        for idx2, a2 in enumerate(a2_values):
            for idx3, b in enumerate(b_values):
                for cost1, prob_cost1 in cost1_dist.items():
                    prob_tmp = 0 
                    for epp, prob_epp in eprimeprime.items():
                        thres_tmp1 = ( b - epp - a1 * cost1 ) / a2
                        thres_tmp2 = ( b + epp + a1 * cost1 ) / a2
                        pfail = pcosthetaSurvival( thres_tmp1 ) + pcosthetaSurvival( thres_tmp2 )
                        prob_tmp += prob_epp * pfail
                    V[idx1, idx2, idx3] += prob_cost1 * prob_tmp
                if V[idx1, idx2, idx3] < 2**-300:  
                    V[idx1, idx2, idx3:] = 0
                    break

    # make interpolation function
    fn = RegularGridInterpolator((a1_values,a2_values,b_values), V)

    # make a list of ciphertext probability and failure probability
    sinthetaSE = sin(thetaSE)
    costhetaSE = cos(thetaSE)

    ########################################################################
    # calculate the failure probability of all ciphertexts in the normc1, normc2, theta1, theta2 list
    ########################################################################

    # for efficiency reasons, prepare computations of P[norms1, norms2]
    s_list = []
    p_list = []
    for norms1, norms2 in itertools.product(dist_norms1, dist_norms2):
        s_list.append([norms1, norms2, 1])
        p_list.append(dist_norms1[norms1] * dist_norms2[norms2])
    s_list = np.squeeze( np.array([s_list]) )
    p_list = np.array([p_list])

    # for efficiency reasons, prepare matrix used to convert [normc1, normc2, theta1, theta2] into [a1, a2, b]
    transformation = np.zeros( (3,3))
    transformation[2,2] = thres

    # loop over all ciphertexts and add failure prob to list
    thelist = []
    for i, tmp in tqdm.tqdm(enumerate(itertools.product(xrange1, xrange2, yrange, yrange)), total=len(xrange1) * len(xrange2) * len(yrange)**2, desc='make list: '):
        normc1, normc2, theta1, theta2 = tmp
        if p_cipher[i] > ACCURACY:
            # finalize matrix used to convert [normc1, normc2, theta1, theta2] into [a1, a2, b]
            transformation[0,0] = normc1 * sinthetaSE * sin(theta1)
            transformation[1,1] = normc2 * sinthetaSE * sin(theta2)
            transformation[0,2] = - normc1 * costhetaSE * cos(theta1)
            transformation[1,2] = - normc2 * costhetaSE * cos(theta2)

            # calculate average failure probability over different norms1, norms2
            tmp = np.dot(s_list, transformation)
            pfail = fn( tmp )
            pfail = np.sum(pfail * p_list)
            
            # add to list
            thelist.append((p_cipher[i], pfail))

    ########################################################################
    # convert the list to alpha beta values
    ########################################################################
    alpha, beta = listToAlphabeta(thelist)

    return alpha, beta