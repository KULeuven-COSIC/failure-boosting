import pickle
import os.path
import numpy as np
from math import sqrt

from gaussian import firstgaussian, moregaussian
from geometric import firstgeometric, firstgeometricsimple, moregeometric
from geometric_uneven import firstgeometric_uneven, moregeometric_uneven

Grover = lambda alpha, MaxDepth: (alpha**0.5) * max(1, alpha**0.5/MaxDepth)

f_alpha = [u'work to generate one weak sample (1/α)', lambda a, b: a]
f_sqrtalpha = [u'work to generate one weak sample (1/√α)', lambda a,b: sqrt(a)]
f_beta = [u'weak ciphertext failure rate (β)', lambda a,b: b]
f_alphabeta = [u'total work to generate a failure (1/αβ)', lambda a, b: a * b**-1]
f_sqrtalphabeta = [u'total work to generate a failure (1/β√α)', lambda a, b: sqrt(a) * b**-1]

f_sqrtalphabetagrover = [u'total work to generate a failure (1/β√α)', lambda a, b, maxdepth: Grover(a, maxdepth) * b**-1]
f_sqrtalpha = [u'work to generate one weak sample (1/√α)', lambda a, b, maxdepth: Grover(a, maxdepth)]

def failureboosting(scheme, method, n_failure, recalc=False, targets=1, beta0=None):
    n = scheme['n'] # dimension of A; We assume square
    n2 = scheme['n2'] # total coefficients in sAs'
    thres = scheme['thres'] # threshold for failures; typically q/4.

    s = scheme['s'] # distribution of s
    sprime = scheme['sprime'] # distribution of s'
    e = scheme['e'] # distribution of e
    eprime = scheme['eprime'] # distribution of e'
    eprimeprime = scheme['eprimeprime'] # distribution of e''

    # load if already generated
    extraname = ""
    if targets > 1 or beta0 != None:
        extraname = f"-T{targets}-beta{int(round(np.log2(beta0)))}"
    name = 'intermediates_fb/'+scheme['name'] + '-' + method + '-' + str(n_failure) + extraname
    if os.path.exists(name + '-' + 'beta.txt') and not recalc:
        alpha = np.loadtxt(name + '-' + 'alpha.txt')
        beta = np.loadtxt(name + '-' + 'beta.txt')
    else:
        if method == 'gaussian':
            if n_failure == 1:
                alpha, beta = firstgaussian(n, n2, thres, s, sprime, e, eprime, eprimeprime)
            else:
                raise Exception("Did not give accurate results because distribution is too skewed and Gaussian assumption does not hold. Code is still available to help future research.") 
                alpha, beta = moregaussian(n, n2, thres, s, sprime, e, eprime, eprimeprime)
        elif method == 'geometric':
            if n_failure == 1:
                alpha, beta = firstgeometric(n, n2, thres, s, sprime, e, eprime, eprimeprime)
            else:
                alpha, beta = moregeometric(n, n2, thres, s, sprime, e, eprime, eprimeprime, n_failure, targets=targets, beta0=beta0)
        elif method == 'geometric-uneven':
            if n_failure == 1:
                alpha, beta = firstgeometric_uneven(n, n2, thres, s, sprime, e, eprime, eprimeprime)
            else:
                alpha, beta = moregeometric_uneven(n, n2, thres, s, sprime, e, eprime, eprimeprime, n_failure, targets=targets, beta0=beta0)
        elif method == 'geometric-simple':
            if n_failure == 1:
                alpha, beta = firstgeometricsimple(n, n2, thres, s, sprime, e, eprime, eprimeprime)
        else:
            raise Exception('not a valid method')

        np.savetxt(name + '-' + 'alpha.txt', alpha)
        np.savetxt(name + '-' + 'beta.txt', beta)

    return alpha, beta