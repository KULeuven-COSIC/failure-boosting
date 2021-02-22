from dist import Dist
from scipy.stats import norm
from math import log
import tqdm
import itertools
import numpy as np

def firstgaussian(n, n2, thres, s, sprime, e, eprime, eprimeprime, n3=1):
    """ calculate alpha-beta list for first failure using gaussian approximation """
    
    def Probl(n, s, sprime, e, eprime):
        """ calculate the distribution of l = var(S^T C)_ijk """

        # get the variance of the secrets
        varS = s.var()
        varE = e.var()

        # calculate the distribution of var(S^T C)_ijk
        abss = Dist({})
        abse = Dist({})

        for i in sprime.keys():
            abss[i**2 * varE] = abss.get(i**2 * varE, 0) + sprime[i]
        for i in eprime.keys():
            abse[i**2 * varS] = abse.get(i**2 * varS, 0) + eprime[i]

        abss = (abss.set_simplify_level(2**12))**n
        abse = (abse.set_simplify_level(2**12))**n
        res = (abss + abse).simplify_int()

        return res


    def Failgivenl(variance):
        """ gaussian approximation of the distribution of the error based on the variance """
        return lambda thres: norm.sf(thres, loc=0, scale=variance**0.5)

    # only abs of eprimeprime is important
    absepp = eprimeprime.abs()

    # calculate the distribution of var(S^T C)_ijk
    probl = Probl(n, s, sprime, e, eprime)

    # calculate the distribution of the failure prob for each value of var(S^T C)_ijk
    thelist = []
    for varl in tqdm.tqdm(probl.keys()):
        eppfail = Dist({})
        fail = Failgivenl(varl)
        for i in absepp.keys():
            eppfail[fail(thres - i) + fail(thres + i)] = absepp[i]
        eppfail = (eppfail.set_simplify_level(2**8))**n2
        for i in eppfail.keys():
            thelist.append((eppfail[i], varl, i))

    # avarage over the values of var(S^T C)_ijk
    faildist = Dist({})
    for pepp, l, failprob in tqdm.tqdm(thelist):
        faildist[failprob] = faildist.get(failprob, 0) + probl[l] * pepp

    # final convolution
    faildist = (faildist.set_simplify_level(2**13))**n3

    ########################################################################
    # convert the list to alpha beta values
    ########################################################################
    thelist = sorted(faildist.keys(), reverse=True)

    # calculate alpha and beta
    alpha, beta = [], []
    alphatmp, betatmp = 0, 0

    for i in thelist:
        alphatmp += faildist[i]
        alpha.append(alphatmp)

        betatmp += i * faildist[i]
        beta.append(betatmp / alphatmp)

    # cutoff values of alpha smaller than 2**-256 (not useful)
    beta2 = [b for (b, a) in zip(beta, alpha)]
    alpha2 = [a**-1 for a in alpha]

    return alpha2, beta2

#####################################################################################
# below is old code that is not functional nor correct. But I leave it here as it might help future research. In general, the approach does not work as the gaussian approximation does not hold (the distribution is too skewed)
#####################################################################################

ACC1 = 2**8
ACC2 = 2**8

def getestimator(thres, s, sprime, e, eprime, eprimeprime, n):
    '''
    determine the estimator of s_i given e'_i
    Notation: s_i -> s, e'_i -> sp
    '''
    print('WARNING: this code will not run nor is it correct. Its included to help potential future research. Use at your own risk')

    svalues = sorted(s.keys())
    evalues = sorted(eprime.keys())

    sgivensp = {i: 1. for i in itertools.product(svalues, evalues)}
    spgivenfail = {i: 1. for i in evalues}

    norms = 0
    normsp = {}

    # make the error distribution without s_i, eprime_i
    seprime = law_product(s, eprime)
    esprime = law_product(sprime, e)

    tmp1 = iter_law_convolution(seprime, n - 1)
    tmp2 = iter_law_convolution(esprime, n)
    tmp = law_convolution(tmp1, tmp2)
    tmp = law_convolution(tmp, eprimeprime)
    tmp = tail_probability(tmp)

    # calcalate probabilities
    for epi in evalues:
        spgivenfail[epi] = eprime[epi]
        prob = 0.
        for si in svalues:
            prob += tmp[int(thres + 1 - epi * si)] * s[si]
        spgivenfail[epi] *= prob
        norms += spgivenfail[epi]

    # calcualte the probabilities
    for (si, epi) in itertools.product(svalues, evalues):
        # P[si]
        sgivensp[(si, epi)] *= s[si]
        # P[fail | si]
        sgivensp[(si, epi)] *= tmp[int(thres + 1 - epi * si)]
        # compute norm
        normsp[epi] = normsp.get(epi, 0) + sgivensp[(si, epi)]


    # normalize
    for (si, epi) in itertools.product(svalues, evalues):
        sgivensp[(si, epi)] /= normsp[epi]

    for epi in evalues:
        spgivenfail[epi] /= norms

    # calculate mean, sigma for a given si
    meanvariancegivensp = {}
    for epi in evalues:
        D = {}
        for si in svalues:
            D[si] = sgivensp[(si, epi)]
        mu = distmean(D)
        var = distvariance(D)
        meanvariancegivensp[epi] = (mu, var)

    return spgivenfail, meanvariancegivensp


def moregaussian(n, n2, thres, s, sprime, e, eprime, eprimeprime, n3=1):
    """ calculate alpha-beta list for first failure using gaussian approximation """
    
    def Probl(clist):
        """ calculate the distribution of var(S^T C)_ijk """
        res = {(0, 0): 1.}
        for n, pc, mu, var in tqdm.tqdm(clist):
            tmp = {}
            for i in pc.keys():
                key = (i*mu, i**2*var)
                tmp[key] = tmp.get(key, 0) + pc[i]
            tmp = iter_law_convolution_simplify_double(tmp, n, ACC1, ACC2)
            res = law_convolution_double(tmp, res, ACC1, ACC2)
            # res = simplifyDistribution_double(res, ACC1, ACC2)

        return clean_dist(res)


    def Failgivenl(mu, variance):
        """ gaussian approximation of the distribution of the error based on the variance """
        return lambda thres: norm.sf(thres, loc=mu, scale=sqrt(variance))

    epgivenfail, meanvariancegivenep = getestimator(thres, s, sprime, e, eprime, eprimeprime, n)
    spgivenfail, meanvariancegivensp = getestimator(thres, e, eprime, s, sprime, eprimeprime, n)

    # generate failing ciphertext
    keys = list(spgivenfail.keys())
    prob = [spgivenfail[k] for k in keys]

    clist = []
    cdict = {}
    for i in range(0, n):
        sp = np.random.choice(keys, p=prob)
        cdict[sp] = cdict.get(sp, 0) + 1

    for sp in cdict.keys():
        clist.append( (cdict[sp], sprime, meanvariancegivensp[sp][0], meanvariancegivensp[sp][1]) )

    # generate failing ciphertext
    keys = list(epgivenfail.keys())
    prob = [epgivenfail[k] for k in keys]

    cdict = {}
    for i in range(0, n):
        ep = np.random.choice(keys, p=prob)
        cdict[ep] = cdict.get(ep, 0) + 1

    for ep in cdict.keys():
        clist.append( (cdict[ep], eprime, meanvariancegivenep[ep][0], meanvariancegivenep[ep][1]) )

    # # REMOVE
    # clist = []
    # clist.append( (n, eprime, 0., distvariance(s) ) )
    # clist.append( (n, sprime, 0., distvariance(e) ) )


    # only abs of eprimeprime is important
    absepp = law_abs(eprimeprime)

    # calculate the distribution of var(S^T C)_ijk
    probl = Probl(clist)

    # REMOVE
    xas = sorted(probl.keys(), key= lambda x: x[1])
    yas = [probl[x] for x in xas]
    xas = [x[1] for x in xas]
    import matplotlib.pyplot as plt
    plt.clf()
    plt.semilogy(xas, yas)
    plt.savefig('multi.pdf')
    plt.show()

    # calculate the distribution of the failure prob for each value of var(S^T C)_ijk
    thelist = []
    for mul, varl in tqdm.tqdm(probl.keys()):
        eppfail = {}
        fail = Failgivenl(mul, varl)
        for i in absepp.keys():
            eppfail[fail(thres - i) + fail(thres + i)] = absepp[i]
        eppfail = iter_law_convolution_simplify(eppfail, n2, 2**8)
        for i in eppfail.keys():
            thelist.append((eppfail[i], (mul, varl), i))

    # avarage over the values of var(S^T C)_ijk
    faildist = {}
    for pepp, key, failprob in tqdm.tqdm(thelist, desc='make list: '):
        faildist[failprob] = faildist.get(failprob, 0) + probl[key] * pepp

    # final convolution
    faildist = iter_law_convolution_simplify(faildist, n3, 2**13)

    # sort the list by failure probability
    thelist = sorted(faildist.keys(), reverse=True)

    # calculate alpha and beta
    alpha, beta = [], []
    alphatmp, betatmp = 0, 0

    for i in thelist:
        alphatmp += faildist[i]
        alpha.append(alphatmp)

        betatmp += i * faildist[i]
        beta.append(betatmp / alphatmp)

    # cutoff values of alpha smaller than 2**-256 (not useful)
    beta2 = [b for (b, a) in zip(beta, alpha)]
    alpha2 = [a**-1 for a in alpha]

    return alpha2, beta2