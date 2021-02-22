import numpy as np
from math import log, pi, cos, sin, sqrt, acos
import tqdm
import itertools
from dist import Dist

# gridpoint that determine accuracy / speed tradeoff
ACCURACY = 2**-300
GRIDPOINTS = 200
GRIDPOINTSS = 100
GRIDPOINTSEPP = 40

def Ptheta(theta, n):
    """Compute probability of a uniformly distributed angle. Not normalized!
    :param theta: angle (radians)
    :param n: dimension"""
    return sin(theta)**(n - 2)  # not normalized!


def PthetaMax(n, m):
    if m == 1:
        return lambda theta: Ptheta(theta, n)

    points = 10000
    yrange = np.linspace(0, float(pi), points)

    f = []
    F = []

    # generate f (probability function)
    # generate F (1 - cumulative function)
    sum = 0
    for theta in yrange:
        prob = Ptheta(theta, n)
        f.append(prob)
        F.append(sum)
        sum += prob

    f = [i / sum for i in f]
    F = [i / sum for i in F]

    fmax = []
    for i, theta in enumerate(yrange):
        fmax.append(m * f[i] * abs(F[i] - 1)**(m - 1))

    fmax = np.array([float(i) for i in fmax])

    return lambda angle: np.interp(float(angle), yrange, fmax)

# calculate survival distribution of cos(theta)

def PcosthetaSurvival(n):
    """Compute the survival function distribution of cos(theta). internal value POINTS can be changes for more accuracy
    :param n: dimension"""
    POINTS = 10000
    thetarange = np.linspace(-1, 1, POINTS)

    pthetasurvival = []
    prob = 0
    for costheta in reversed(thetarange):
        prob += float(sin(acos(costheta))**(n-2))
        pthetasurvival.append(float(prob))

    pthetasurvival = [i / prob for i in reversed(pthetasurvival)]

    pthetasurvival = np.array(pthetasurvival)
    return lambda angle: np.interp(float(angle), thetarange, pthetasurvival)

def listToAlphabeta(thelist):
    """convert a list of ciphertexts with a given probability and a given failure probability into alpha and beta values"""
    # cleanup by combining similar items
    thelist.sort(key=lambda x: x[0])

    dict_list = {}
    for Pcipher, Pfail in thelist:
        key = float(np.log2(Pfail))
        key = round(key, 2)
        if key in dict_list:
            Pcipher2, Pfail2 = dict_list[key]
            Pciphertot = Pcipher + Pcipher2
            if Pciphertot < ACCURACY:
                dict_list[key] = ( Pciphertot, Pfail )
            else:
                dict_list[key] = (Pcipher + Pcipher2, (Pfail * Pcipher + Pfail2 * Pcipher2) / (Pcipher + Pcipher2) )
        else:
            dict_list[key] = (Pcipher, Pfail)
    
    thelist= list(dict_list.values())


    # sort the list in descending failure probability
    thelist.sort(key=lambda x: x[1], reverse=True)

    # calculate alpha and beta
    alpha2, beta2 = [], []
    alphatmp, betatmp = 0, 0

    for i in thelist:
        alphatmp += i[0]
        alpha2.append(alphatmp)

        betatmp += i[0] * i[1]
        beta2.append(betatmp / alphatmp)

    # get into same form as plots from decryption failure paper
    beta = [float(b) for (b, a) in zip(beta2, alpha2)]
    alpha = [float(a**-1) for a in alpha2]

    return alpha, beta


def calcFailgivennormctNoInfo(normc, normsdist, epp, thres, n2, pcosthetaSurvival):
    """ Calculate the failure probability given a certain norm of C.
    Takes into account distribution norm of S and e'' distribution.
    Approximations: geometric approximation
    :param normc: norm of C
    :param normsdist: distribution of norm of S
    :param epp: distribution of e''
    :param thres: failure threshold (typically q/4)
    :param n2: coefficients in transmitted messages
    :param pcosthetaSurvival: survival function of theta
    """
    prob = 0
    for norms in normsdist.keys():
        probtmp = 0
        for g in epp.keys():
            thres1 = (thres - g) / normc / norms
            thres2 = (thres + g) / normc / norms

            probtmp += epp[g] * (pcosthetaSurvival(thres1) +
                                 pcosthetaSurvival(thres2))

        # take into account n2 coefficients
        if probtmp > 2**-10:
            probtmp = 1 - (1 - probtmp)**n2
        else:
            probtmp = probtmp * n2

        prob += probtmp * normsdist[norms]

    return prob


def firstgeometric(n, n2, thres, s, sprime, e, eprime, eprimeprime):
    """ calculate alpha-beta list for first failure using geometric approximation. 
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
    # Make a grid of possible normc1, normc2 combinations and calculate the probability
    ########################################################################
    normcdist = sprime.norm(n, eprime).simplify(GRIDPOINTS).normalize()
    xrange = normcdist.keys()

    # get the probability of ciphertexts over the grid
    Pcipher = []
    sumprob = 0
    for normc in tqdm.tqdm(xrange, total=len(xrange), desc='P[ciphertext]: '):
        prob = normcdist[normc]
        Pcipher.append(prob)
        sumprob += prob

    # normalize the probability
    Pcipher = [i / sumprob for i in Pcipher]

    ########################################################################
    # calculate the failure probability of all ciphertexts in the normc1, normc2 list
    ########################################################################

    # calculate norms
    normsdist = s.norm(n, e).simplify(GRIDPOINTSS).normalize()
    eprimeprime = eprimeprime.simplify(GRIDPOINTSEPP)

    # make a list of ciphertext probability and failure probability
    pcosthetaSurvival = PcosthetaSurvival(2*n)
    thelist = []
    for i, normc in tqdm.tqdm(enumerate(xrange), total=len(xrange), desc='make list: '):
        if Pcipher[i] > ACCURACY:
            fail = calcFailgivennormctNoInfo(
                normc, normsdist, eprimeprime, thres, n2, pcosthetaSurvival)
            thelist.append((Pcipher[i], fail))

    ########################################################################
    # convert the list to alpha beta values
    ########################################################################
    alpha, beta = listToAlphabeta(thelist)

    return alpha, beta


def calcFailgivennormctNoInfoSimple(normc, norms, thres, n2, pcosthetaSurvival):
    """ Calculate the failure probability given a certain norm of C and norm of S.
    Approximations: 
    - geometric approximation
    - does not take into account e'' distribution
    :param normc: norm of C
    :param norms: norm of S
    :param thres: failure threshold (typically q/4)
    :param n2: coefficients in transmitted message
    :param pcosthetaSurvival: survival function of theta
    """
    thres = (thres / normc / norms)

    # deal with out of scope values
    prob = 2 * pcosthetaSurvival(thres)

    if prob > 2**-10:
        prob = 1 - (1 - prob)**n2
    else:
        prob = prob * n2

    return prob


def firstgeometricsimple(n, n2, thres, s, sprime, e, eprime, eprimeprime):
    """ calculate approximate alpha-beta list for first failure using geometric approximation. This is essentially EC:DanRosVir20 method
    Approximations: 
    - geometric approximation, 
    - equal distribution of s and e and equal distribution of s' and e', 
    - approximates some distribution to GRIDPOINTS values for efficiency reasons, 
    - does not take into account e'' distribution
    - does not take into account S distribution 
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
    # Make a grid of possible normc1, normc2 combinations and calculate the probability
    ########################################################################
    normcdist = sprime.norm(n, eprime).simplify(GRIDPOINTS).normalize()
    xrange = normcdist.keys()

    # get the probability of ciphertexts over the grid
    Pcipher = []
    sumprob = 0
    for normc in tqdm.tqdm(xrange, total=len(xrange), desc='P[ciphertext]: '):
        prob = normcdist[normc]
        Pcipher.append(prob)
        sumprob += prob

    # normalize the probability
    Pcipher = [i / sumprob for i in Pcipher]

    ########################################################################
    # calculate the failure probability of all ciphertexts in the normc1, normc2 list
    ########################################################################
    normsdist = s.norm(n, e).simplify(GRIDPOINTSS).normalize()
    norms = normsdist.mean()

    # make a list of ciphertext probability and failure probability
    pcosthetaSurvival = PcosthetaSurvival(2*n)
    thelist = []
    for i, normc in tqdm.tqdm(enumerate(xrange), total=len(xrange), desc='make list: '):
        if Pcipher[i] > ACCURACY:
            fail = calcFailgivennormctNoInfoSimple(
                normc, norms, thres, n2, pcosthetaSurvival)
            thelist.append((Pcipher[i], fail))

    ########################################################################
    # convert the list to alpha beta values
    ########################################################################
    alpha, beta = listToAlphabeta(thelist)

    return alpha, beta


def calcErrorgivennormct(normc, normsdist, epp, thres, pcosthetaSurvival, thetaCE, thetaSE):
    """ Calculate the !error! probability given a certain norm of C and angle with E (directional failure boosting).
    Approximations: 
    - geometric approximation
    :param normc: norm of C
    :param norms: distribution of norm of S
    :param epp: distribution of e''
    :param thres: failure threshold (typically q/4)
    :param pcosthetaSurvival: survival function of theta
    :param thetaCE: theta_CE
    :param thetaSE: theta_SE
    """
    prob = 0
    for norms in normsdist.keys():
        probtmp = 0
        for g in epp.keys():
            # filter out the zero case
            if (normc * norms * sin(thetaCE) * sin(thetaSE)) == 0:
                if abs(normc * norms * cos(thetaCE) * cos(thetaSE) + g) > thres:
                    probtmp += epp[g]
                    continue
                else:
                    continue

            # calculate thresholds
            thres1 = (thres - g - normc * norms * cos(thetaCE) *
                      cos(thetaSE)) / (normc * norms * sin(thetaCE) * sin(thetaSE))
            thres2 = (thres + g + normc * norms * cos(thetaCE) *
                      cos(thetaSE)) / (normc * norms * sin(thetaCE) * sin(thetaSE))

            probtmp += epp[g] * (pcosthetaSurvival(thres1) +
                                 pcosthetaSurvival(thres2))

        prob += probtmp * normsdist[norms]

    return prob


def normSgivenft(normsdist, normcdistft, epp, n2, thres, targets, pcosthetaSurvival):
    # determine P[F | ft] and P[F | norms, minnormc0]
    PFgivenSft = Dist({})
    for norms in normsdist.keys():
        PFgivenSft[norms] = 0
        for normc in normcdistft.keys():
            PFgivenSft[norms] += normcdistft[normc] * calcFailgivennormctNoInfo(
                normc, {norms: 1.}, epp, thres, n2, pcosthetaSurvival)

    PFgivenSft = PFgivenSft.normalize()
    PFgivenft = sum([PFgivenSft[norms]*normsdist[norms]
                     for norms in normsdist])

    # calculate new normsdist
    normsdistNew = {}
    for norms in normsdist.keys():
        normsdistNew[norms] = normsdist[norms] * PFgivenSft[norms] / \
            (PFgivenSft[norms] + (targets-1) * PFgivenft)

    return Dist(normsdistNew).normalize()


def moregeometric(n, n2, thres, s, sprime, e, eprime, eprimeprime, n_failure, targets=1, beta0=None):
    """calculate directional failure boosting using geometric approach"""

    ########################################################################
    # Make a grid of possible normc1, normc2, theta_CE1, theta_CE2 combinations
    ########################################################################
    normcdist_full = sprime.norm(n, eprime)
    normcdist = normcdist_full.simplify(GRIDPOINTS).normalize()
    xrange = normcdist.keys()

    pthetamax = PthetaMax(2 * n, n2)
    angledist = Dist({ angle : pthetamax(angle) for angle in np.linspace(0, float(pi), 10000) })
    angledist = angledist.normalize().simplify(GRIDPOINTS)
    yrange = angledist.keys()

    ########################################################################
    # calculate the probability of all ciphertexts in the grid
    ########################################################################
    Pcipher = []
    sumprob = 0
    for normc, theta in tqdm.tqdm(itertools.product(xrange, yrange), total=len(xrange) * len(yrange), desc='P[ciphertext]: '):
        prob = angledist[theta] * normcdist[normc]
        Pcipher.append(prob)
        sumprob += prob

    # normalize the probability
    Pcipher = [i / sumprob for i in Pcipher]

    # calculate angle survival function
    pcosthetaSurvival = PcosthetaSurvival(2*n)

    ########################################################################
    # Determine the mean norm of c1, c0 in first stage
    ########################################################################
    if beta0 == None:
        minnormc0 = 0
    else:
        normcdist = sprime.norm(n, eprime)
        distance = 2**300
        tailprobnormc = normcdist.sf()
        for normc in normcdist.keys():
            prob = tailprobnormc(normc)
            if np.abs(prob - beta0) < distance:
                distance = np.abs(prob - beta0)
                minnormc0 = normc

    normcdistft = Dist(
        {normc: normcdist_full[normc] for normc in normcdist_full if normc >= minnormc0})
    normcdistft = normcdistft.normalize().simplify(GRIDPOINTS)

    ########################################################################
    # Take into account weak secret if multiple targets are attacked in first stage
    ########################################################################
    if targets > 1:
        normsdist_apriori = s.norm(n, e)
        normsdist = normSgivenft(
            normsdist_apriori, normcdistft, eprimeprime, n2, thres, targets, pcosthetaSurvival)
    else:
        normsdist = s.norm(n, e)

    ########################################################################
    # Determine theta_SE
    ########################################################################
    meanNorms = normsdist.mean()
    meanNormc = normcdistft.mean()

    thetaCS = acos(thres / meanNorms / meanNormc)
    costhetaSE = cos(thetaCS) / sqrt(cos(thetaCS)**2 +
                                     sin(thetaCS)**2 / (n_failure-1))
    thetaSE = float(acos(costhetaSE))
    pcosthetaSurvival = PcosthetaSurvival(2*n - 1)

    # simplify normsdist for efficiency
    normsdist = normsdist.simplify(GRIDPOINTSS)
    eprimeprime = eprimeprime.simplify(GRIDPOINTSEPP)

    ########################################################################
    # calculate error probability if |cos(theta_j)| < |cos(theta_i)| 
    # this is the residual error probability of ciphertexts at different locations
    ########################################################################
    pcosthetaSurvival = PcosthetaSurvival(2*n - 1)
    failOther = {}
    for normc in tqdm.tqdm(xrange, desc='P[F | |cos(theta_j)| < |cos(theta_i)|]: '):
        sumthetaprob = 0
        sumfailprob = 0
        for theta in sorted(yrange, key=lambda x: abs(x - pi / 2)):
            ptheta = Ptheta(theta, 2 * n)
            if ptheta < 2**-300:
                failOther[(normc, theta)] = float(sumfailprob / sumthetaprob)
            else:
                sumthetaprob += ptheta
                sumfailprob += ptheta * \
                    calcErrorgivennormct(
                        normc, normsdist, eprimeprime, thres, pcosthetaSurvival, theta, thetaSE)
                if sumthetaprob > 0:
                    failOther[(normc, theta)] = float(
                        sumfailprob / sumthetaprob)
                else:
                    failOther[(normc, theta)] = 0

    ########################################################################
    # calculate the probability of all ciphertexts in the grid
    ########################################################################
    thelist = []
    for i, tmp in tqdm.tqdm(enumerate(itertools.product(xrange, yrange)), total=len(xrange) * len(yrange), desc='make list: '):
        normc, theta = tmp
        if Pcipher[i] > ACCURACY:
            Pfailmaxcosine = calcErrorgivennormct(
                normc, normsdist, eprimeprime, thres, pcosthetaSurvival, theta, thetaSE)
            Pfailothercosine = failOther[(normc, theta)]
            if Pfailothercosine > 2**-10:
                prob = 1 - (1 - Pfailmaxcosine) * \
                    (1 - Pfailothercosine)**(n2-1)
            else:
                prob = Pfailmaxcosine + \
                    (n2-1) * Pfailothercosine - \
                    Pfailmaxcosine * (n2-1) * Pfailothercosine
            thelist.append((Pcipher[i], prob))

    ########################################################################
    # convert the list to alpha beta values
    ########################################################################
    thelist.sort(key=lambda x: x[1], reverse=True)

    alpha, beta = listToAlphabeta(thelist)

    return alpha, beta
