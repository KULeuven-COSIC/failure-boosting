import numpy as np
from math import log
import matplotlib.pyplot as plt
import copy
import tqdm
import os

# from plot import plot
from failureboosting import failureboosting, f_alpha, f_sqrtalpha, f_beta, f_alphabeta, f_sqrtalphabeta, Grover




def firstfailureonly(scheme, maxqueries, targets, recalc=False, method='geometric-uneven', limitct = float('inf'), maxdepth = float('inf')):
    alpha1, beta1 = failureboosting(scheme=scheme, method = "gaussian", n_failure = 1, recalc=recalc)
    querylimit = maxqueries * targets
    sqrtalphabeta1 = [ (Grover(a, maxdepth)/b, b) for (a, b) in zip(alpha1, beta1) if 1/b < querylimit and a/b/targets < limitct ]
    if len(sqrtalphabeta1) == 0:
        return float('inf'), float('inf')
    w1, q1 = min(sqrtalphabeta1, key = lambda x: x[0])
    return w1, q1


def multitarget(scheme, maxqueries, targets, recalc=False, simpleMultitarget=False, addFail3 = True, plotit = True, legendloc='upper right', method='geometric-uneven', limitct = float('inf'), maxdepth = float('inf')):

    ## find first failure using regular failure boosting
    alpha1, beta1 = failureboosting(scheme, method = "gaussian", n_failure = 1, recalc=False)

    # divide the maximal number of queries over the first three failures
    maxqueries = maxqueries / 3

    # find normc in normal failure boosting for reference (ideally this should be adapted in the loop, but it is too costly)
    querylimit1 = maxqueries * targets
    # print(np.log2(querylimit1))
    sqrtalphabeta1 = [ (Grover(a, maxdepth)/b, b) for (a, b) in zip(alpha1, beta1) if 1/b < querylimit1 ]
    if len(sqrtalphabeta1) == 0:
        return float('inf'), float('inf'), float('inf'), 0, 0, 0, {'overshoot': 0, 'paramsSearch2':{}}
    w1, q1 = min(sqrtalphabeta1, key = lambda x: x[0])

    ## find second failure optimally
    alpha2, beta2 = failureboosting(scheme, method = method, n_failure = 2, recalc=recalc, targets=targets, beta0=q1)

    paramsSearch2 = {'targets' : targets, 'beta0' : q1 }
    
    work1 = float('inf')
    work2 = float('inf')
    queries1 = 0
    queries2 = 0
    overshootoptimal = 0

    # prepare a list with overshoot values
    if simpleMultitarget:
        overshootlist = [0] # no overshoot for simple multitarget
    else:
        maxovershoot = int( log(targets, 2) )
        overshootlist = np.linspace(0, maxovershoot, 200) # between 0 and targets overshoot for levelled multitarget

    # loop over all overshoot values and find optimal value of overshoot
    for overshoot in tqdm.tqdm( (int(2**i) for i in overshootlist), total=len(overshootlist), leave=None):
        # search for overshoot number of failures
        querylimit1 = maxqueries * targets / overshoot
        sqrtalphabeta1 = [ (Grover(a, maxdepth)/b, b) for (a, b) in zip(alpha1, beta1) if 1/b < querylimit1 and a/b/targets*overshoot < limitct ]

        if len(sqrtalphabeta1) == 0:
            continue
        w1, q1 = min(sqrtalphabeta1, key= lambda x: x[0])
        w1 = w1 * overshoot
        q1 = q1 / overshoot

        # find second failure given overshoot number of targets
        querylimit2 = maxqueries * overshoot
        sqrtalphabeta2 = [ (Grover(a, maxdepth)/b, b) for (a, b) in zip(alpha2, beta2) if 1/b < querylimit2 and a/b/overshoot < limitct ]
        if len(sqrtalphabeta2) == 0:
            continue
        w2, q2 = min(sqrtalphabeta2, key= lambda x: x[0])

        if w1 + w2 < work1 + work2:
            work1, queries1 = w1, q1
            work2, queries2 = w2, q2
            overshootoptimal = overshoot

    # ## find third failure optimally
    if addFail3:
        alpha3, beta3 = failureboosting(scheme, method = method, n_failure = 3, recalc=recalc, targets=targets, beta0=q1)

        # find optimal failure boosting strategy
        querylimit = maxqueries

        sqrtalphabeta3 = [ (a**0.5/b, b) for (a, b) in zip(alpha3, beta3) if 1/b < querylimit and a/b < limitct ]
        if len(sqrtalphabeta3) == 0:
                return 2**300
        work3, queries3 = min(sqrtalphabeta3, key= lambda x: x[0])
    else:
        work3 = 0
        queries3 = float("inf")

    return work1, work2, work3, queries1, queries2, queries3, {'overshoot': overshootoptimal, 'paramsSearch2':paramsSearch2}