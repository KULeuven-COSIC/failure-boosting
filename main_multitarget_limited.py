import numpy as np
from failureboosting import failureboosting, f_alpha, f_sqrtalpha, f_beta, f_alphabeta, f_sqrtalphabeta, f_sqrtalphabetagrover
from multitargetattack import multitarget, firstfailureonly
from dist import Dist
from math import log, isinf

log2i = lambda x : int(log(x,2)) if not isinf(log(x,2)) else log(x,2)

import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14}) # for small paper figures
plt.rcParams.update({'font.size': 12})

from NISTschemes import LightSaber, Saber, FireSaber, Kyber512, Kyber768, Kyber1024, uLightSaber, uSaber, uFireSaber 
from NISTschemes import LightSaberTplus1, SaberTplus1, SaberTplus2, FireSaberTplus1, FireSaberTplus2, FireSaberTplus3

def multitarget_limited(scheme, maxqueries, targets, simpleMultitarget=True, limitct = float('inf'), maxdepth = float('inf'), recalc=False, addFail3=False, method='geometric-uneven', legendloc = 'upper left', outname = ''):
    ###############################
    ## print all the information ##
    ###############################
    print('------------------------------------------')
    print(f'\033[1m{scheme["name"]}')
    print(f'maxqueries: 2^{log2i(maxqueries)}, targets: 2^{log2i(targets)}, method: {method}')
    attack = 'simple' if simpleMultitarget else 'levelled'
    print(f'{attack}, limitct: 2^{log2i(limitct)}, maxdepth: 2^{log2i(maxdepth)} \033[0m')


    # failure probability
    print('')
    dist = ( scheme['e'] * scheme['sprime'] + scheme['eprime'] * scheme['s'] )**scheme['n'] + scheme['eprimeprime']
    pfail = dist > scheme['thres']
    print(f'failure probability: 2^{np.log2(pfail*256)}')

    # only first failure attack
    print('')
    w1, q1 = firstfailureonly(scheme=scheme, maxqueries=maxqueries, targets=targets, recalc=recalc, method=method, limitct = limitct, maxdepth = maxdepth)
    print(f'only first failure cost: {np.log2(w1)}, queries: {np.log2(q1)}')

    # full attack cost
    print('')
    work1, work2, work3, queries1, queries2, queries3, other = multitarget(scheme=scheme, maxqueries=maxqueries, targets=targets, simpleMultitarget=simpleMultitarget, recalc=recalc, addFail3=addFail3, method=method, legendloc = 'upper left', limitct =limitct, maxdepth=maxdepth)
    overshootoptimal = other['overshoot']
    paramsSearch2 = other['paramsSearch2']

    # determine whether attack was possible
    if queries1 == 0:
        AttackPossible = False
    else: 
        AttackPossible = True

    if not AttackPossible:
        print('attack not possible')
    else: 
        print(f'first failure cost: {np.log2(work1)}, queries: {np.log2(queries1)}')
        print(f'second failure cost: {np.log2(work2)}, queries: {np.log2(queries2)}')

        print(f"total optimal work: {np.log2(work1 + work2 + work3)}, {np.log2([work1, work2, work3])}")
        print(f"total queries: {np.log2(np.float64(1)/queries1 + np.float64(1)/queries2 + np.float64(1)/queries3)}, {np.log2([queries1, queries2, queries3])}")

    ##############################
    ## plot all the information ##
    ##############################
    scattersize = 70

    colors = ['tab:blue', 'tab:orange', 'tab:green']

    plt.xlabel(f_beta[0])
    plt.ylabel(f_sqrtalphabeta[0])

    # plot the attack curves
    toplot = [(scheme, 'gaussian', 1, f"failure: 1", {}, {'color' : colors[0]} ), 
                (scheme, method, 2, f"failure: 2", paramsSearch2, {'color' : colors[1]} ), 
                ]
    if addFail3:
        toplot.append( (scheme, method, 3, f"failure: 3", paramsSearch2, {'color' : colors[2]} ) ) 

    # print the curves   
    for scheme, method, n_failure, label, fboptions, plotoptions in toplot:
        alpha, beta = failureboosting(scheme, method, n_failure, **fboptions)

        x = np.array([ f_beta[1](a,b) for a,b in zip(alpha, beta) ])
        y = np.array([ f_sqrtalphabeta[1](a,b) for a,b in zip(alpha, beta) ])
        plt.loglog(x, y, label=label, basex=2, basey=2, **plotoptions)

    # add changes of curves due to Grover depth limit
    if maxdepth < float('inf'):
        for scheme, method, n_failure, label, fboptions, plotoptions in toplot:
            alpha, beta = failureboosting(scheme, method, n_failure, **fboptions)

            x = np.array([ f_beta[1](a,b) for a,b in zip(alpha, beta) ])
            y = np.array([ f_sqrtalphabetagrover[1](a,b, maxdepth) for a,b in zip(alpha, beta) ])
            plt.loglog(x, y, basex=2, basey=2, linestyle = 'dotted', **plotoptions)

    # plot the attack points if possible
    if AttackPossible:
        if work1 < 2**300:
            plt.scatter( queries1 * overshootoptimal, work1 / overshootoptimal, marker = 'X', s=scattersize)
            plt.scatter( queries2, work2, marker = 'X', s=scattersize)
            if addFail3:
                plt.scatter(queries3, work3, marker = 'X', s=scattersize)

            if not simpleMultitarget:
                plt.scatter( queries1, work1, marker = 'o', s=scattersize, color='tab:blue')

    # plot the maxquery constraints
    plt.axvline(x=1/maxqueries, color='tab:grey', linestyle='dotted')
    plt.axvline(x=1/maxqueries/targets, color='tab:grey', linestyle='dotted')

    # plot limitct contraint
    if limitct < float('inf'):
        xleft, xright, ybottom, ytop = plt.axis()
        betavalues = 2**np.linspace(np.log2(xleft), np.log2(xright), 1000)
        xvalues = []
        yvalues = []
        for x in betavalues:
            y = (limitct / x * max(1, x**-1 / maxqueries) )**0.5
            if y < ytop and y > ybottom:
                xvalues.append(x)
                yvalues.append(y)
        plt.loglog(xvalues, yvalues, color='tab:red', linestyle='dashdot', basex=2, basey=2, label=f'(αβ)⁻¹ = 2^{int(log(limitct,2))} T^¹')

    # # JUST TO HAVE NICE PLOTS IN PAPER
    # ax = plt.gca()
    # ax.set_ylim([2.**102,2.**124])
    # ax.set_xlim([2.**-120,2.**-54])

    plt.title(f'{scheme["name"]}: {attack}\nqlimit: 2^{log2i(maxqueries)}, T^0: 2^{log2i(targets)}, |M|: 2^{log2i(limitct)}, Dmax: 2^{log2i(maxdepth)}')
    plt.legend(loc=legendloc)
    plt.savefig(outname + '.png')
    # plt.show()
    plt.clf()
    plt.close('all')


def main():
    # equivalent AES levels
    LightSaber['AES'] = 2**128
    Saber['AES'] = 2**192
    FireSaber['AES'] = 2**256
    Kyber512['AES'] = 2**128
    Kyber768['AES'] = 2**192
    Kyber1024['AES'] = 2**256
    uLightSaber['AES'] = 2**128
    uSaber['AES'] = 2**192
    uFireSaber['AES'] = 2**256

    SaberTplus1['AES'] = 2**192
    SaberTplus2['AES'] = 2**192
    FireSaberTplus1['AES'] = 2**256
    FireSaberTplus2['AES'] = 2**256
    FireSaberTplus3['AES'] = 2**256

    Schemes = []
    Schemes = Schemes + [LightSaber, Saber, FireSaber]
    Schemes = Schemes + [Kyber512, Kyber768, Kyber1024]
    Schemes = Schemes + [SaberTplus1, SaberTplus2, FireSaberTplus1, FireSaberTplus2, FireSaberTplus3]

    Schemes= [Saber]

    method = 'geometric-uneven'
    maxqueries = 2**64
    targets = 2**64
    maxdepths = [2**96, float('inf')]
    limitcts = [2**256, 'AES', float('inf')]
    simpleMultitargets = [True, False]

    for scheme in Schemes:
        limitcts_tmp = [ i if i != 'AES' else scheme['AES'] for i in limitcts ]
        for maxdepth in maxdepths:
            for limitct in limitcts_tmp:
                for simpleMultitarget in simpleMultitargets:
                    outname = f'fig_mult/MT-{scheme["name"]}-maxqueries{log2i(maxqueries)}-targets{log2i(targets)}-Simple{simpleMultitarget}-limited-{log2i(limitct)}-{log2i(maxdepth)}'
                    multitarget_limited(scheme=scheme, maxqueries=maxqueries, targets=targets, simpleMultitarget=simpleMultitarget, limitct=limitct, maxdepth=maxdepth, recalc=False, addFail3=False, method=method,legendloc='upper left', outname=outname)

if __name__ == '__main__':
    main()