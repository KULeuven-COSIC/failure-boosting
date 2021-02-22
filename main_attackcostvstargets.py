import numpy as np
from failureboosting import failureboosting, f_alpha, f_sqrtalpha, f_beta, f_alphabeta, f_sqrtalphabeta, f_sqrtalphabetagrover
from multitargetattack import multitarget
from dist import Dist
import tqdm
import os
from math import log, isinf
import random

from NISTschemes import LightSaber, Saber, FireSaber, Kyber512, Kyber768, Kyber1024, uLightSaber, uSaber, uFireSaber
from NISTschemes import LightSaberTplus1, SaberTplus1, SaberTplus2, FireSaberTplus1, FireSaberTplus2, FireSaberTplus3

log2i = lambda x : int(log(x,2)) if not isinf(log(x,2)) else log(x,2)

import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
colors = plt.rcParams["axes.prop_cycle"].by_key()['color']
fig, ax = plt.subplots()

def main():
    # setup of plot
    recalc = False # recalc everything (costly!)
    addFail3 = False # do you want to take the cost of finding the third failure into account (this is costly and usually doesn't change the cost calculation)
    offsetvalue = 1.5 # slightly separate lines in plot for visibility, set to 1 to remove

    # attack constraints
    maxqueries = 2**64
    targets = [0,5,10,15,20,25,30,35,40,45,50,55,60,64]
    limitct = 2**256 # |M|, number of possible ciphertexts
    maxdepth = 2**96 # quantum computer depth

    # schemes to plot
    schemes = []
    schemes += [LightSaber, Saber, FireSaber]
    schemes += [Kyber512, Kyber768, Kyber1024]

    # what to plot
    plotSimple = False # simple attack
    plotNoConstraints = True # no limitct or maxdepth constraints
    plotSimpleNoConstraints = False # no limitct or maxdepth constraints + simple attack
    plotMoreConstraints = True # limitct determined by AES level
    plotSimpleMoreConstraints = False # limitct determined by AES level + simple attack    

    # include necessary plots (in order from less to more constraints)
    plots = []

    if plotNoConstraints:
        plots.append({'limitct' : float('infinity'), 'maxdepth' : float('infinity'), 'simpleMultitarget' : False}) 
    if plotSimpleNoConstraints:
        plots.append({'limitct' : float('infinity'), 'maxdepth' : float('infinity'), 'simpleMultitarget' : True}) 

    plots.append({'limitct' : limitct, 'maxdepth' : maxdepth, 'simpleMultitarget' : False})
    if plotSimple:
        plots.append({'limitct' : limitct, 'maxdepth' : maxdepth, 'simpleMultitarget' : True}) 

    if plotMoreConstraints:
        plots.append({'limitct' : 'AES', 'maxdepth' : maxdepth, 'simpleMultitarget' : False}) 
    if plotSimpleMoreConstraints:
        plots.append({'limitct' : 'AES', 'maxdepth' : maxdepth, 'simpleMultitarget' : True}) 

    # loop over all schemes
    for scheme in tqdm.tqdm(schemes):
        # prepare print
        print('------------------------------------------')
        print(f'\033[1m{scheme["name"]}\033[0m')
        print()

        # prepare plot
        plt.rcParams.update({'font.size': 12})
        colors = plt.rcParams["axes.prop_cycle"].by_key()['color']

        # for the no-attack line on the plot
        maxwork = 1
        noattack = []
        
        # for the slight offset for visibility reasons
        offset = offsetvalue ** (-len(plots)/ 2 + 0.5)
        
        # loop over the different constraint scenarios
        for args in plots:
            # unpack
            limitct_tmp, maxdepth, simpleMultitarget = args.values()
            attack = 'simple' if simpleMultitarget else 'levelled'

            # include AES level if different
            if limitct_tmp == 'AES': 
                limitct_tmp = scheme['AES'] 
                if limitct == limitct_tmp:
                    continue
            args_tmp = args.copy()
            args_tmp['limitct'] = limitct_tmp

            # load if already calculated, otherwise calculate
            name = f'intermediates_AvsT/AvsT-{scheme["name"]}-maxqueries{log2i(maxqueries)}-targets{targets}-{simpleMultitarget}-limited-{log2i(limitct_tmp)}-{log2i(maxdepth)}'
            if os.path.exists(f'{name}-work.txt') and not recalc:
                work = np.loadtxt(f'{name}-work.txt')
                queries = np.loadtxt(f'{name}-queries.txt')
            else: 
                work = []
                queries = []
                
                # determine the cost for each number of targets
                for i in tqdm.tqdm(targets):
                    targets_tmp = 2**i
                    w1, w2, w3, q1, q2, q3, _ = multitarget(scheme=scheme, maxqueries=maxqueries, targets=targets_tmp, recalc=recalc, addFail3=addFail3, **args_tmp)
                    work_tmp = w1+w2+w3
                    queries_tmp = float(np.float64(1)/q1 + np.float64(1)/q2 + np.float64(1)/q3)
                    work.append(work_tmp)
                    queries.append(queries_tmp)

                # save calculated results
                np.savetxt(f'{name}-work.txt', work)
                np.savetxt(f'{name}-queries.txt', queries)

            # print
            print(f'maxqueries: 2^{log2i(maxqueries)}, method: {attack}')
            print(f'limitct: 2^{log2i(limitct_tmp)}, maxdepth: 2^{log2i(maxdepth)}')
            print('')
            print(f'targets: {targets}')
            print(f'work: {[log2i(w) for w in work]}')
            print(f'queries: {[log2i(q) for q in queries]}')
            print(f'\n\n')
            
            # plot
            targetspow = [2**i for i in targets]
            color=colors.pop(0)
            label=f'{attack}, |M|: 2^{log2i(limitct_tmp)}, Dmax: 2^{log2i(maxdepth)}'

            plt.loglog(targetspow, np.array(work)*offset, label=label, color=color, basey=2, basex=2)
            plt.loglog(targetspow, np.array(queries)*offset, linestyle='dashed', color=color, basey=2, basex=2)
            
            # save for the noattack plot
            maxwork_tmp =  2.**10 * max([1] + [ i for i in work if i<2**300 ])
            maxwork = max(maxwork, maxwork_tmp)
            work_inf = [offset if i>=2**300 else float('inf') for i in work ]
            noattack.append( (work_inf, color) )
            offset *= offsetvalue

        # no-attack plot can be done now that we know maxwork
        for work_inf, color in noattack:  
            plt.loglog(targetspow, np.array(work_inf)*maxwork, color=color, basey=2, basex=2)
        locs, labels = plt.yticks()
        locs = [ i for i in locs if i < maxwork/2**5 ]
        labels = [ log2i(i) for i in locs ]
        locs.append(maxwork)
        labels.append('no attack')
        plt.yticks(locs, labels)

        # and finalize the plot
        plt.legend()
        # plt.title(scheme['name'])
        plt.xlabel('targets')
        plt.ylabel('attack cost')
        plt.tight_layout()
        plt.savefig(f'fig_AvsT/{scheme["name"]}-maxqueries{log2i(maxqueries)}-targets{targets}-attackcostvstargets.png')
        plt.clf()

if __name__ == '__main__':
    main()