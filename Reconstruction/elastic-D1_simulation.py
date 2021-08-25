#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 14:15:46 2017

@original author: dkita
@implemented: joshualchen
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from sklearn.metrics import r2_score
from sklearn.linear_model import ElasticNet
import os
import sys
import scipy.optimize as so
import timeit
import time
from cvxopt import matrix, solvers
import pickle
solvers.options['show_progress']=False

font = {'weight' : 'normal',
        'size'   : 14}
labelsize = 16
ticksize = 14
matplotlib.rc('font', **font)

""" First import A matrices & y-values
"""

pickle_file = 'stored/803_64_2_[fsrmean_fsrmiddle_T1_T2]_pathlengths_statetable_stateheader.pickle'

# loading the characterization matrix (results_T) and the path_lengths
with open(pickle_file, 'rb') as f:
    [[FSR_mean, FSR_middle, results_T1, results_T2], path_lengths, state_table, state_header] = pickle.load(f)

sort_order = np.argsort(path_lengths)  # getting path lengths order to organize characterization matrix

# developing the characterization matrix
characterization_matrix_1_flipped = results_T1['mode 1 transmission']
# A1 = np.square(np.absolute(np.transpose(characterization_matrix_1_flipped)))[sort_order] * 10 ** 3
A1 = np.square(np.absolute(np.transpose(characterization_matrix_1_flipped))) * 10 ** 3
characterization_matrix_2_flipped = results_T2['mode 1 transmission']
# A2 = np.square(np.absolute(np.transpose(characterization_matrix_2_flipped)))[sort_order] * 10 ** 3
A2 = np.square(np.absolute(np.transpose(characterization_matrix_2_flipped))) * 10 ** 3
# A1 = A2
""" Choose type of signals we want to use:
"""
options = ["MZI1"]#, "MZI2", "MZI1MZI2"]

""" Elastic-D1 reconstruction on each type of signal listed in 'options' above
"""
for signal in options:
    start = time.time()
    """ Uncomment below for broadband reconstruction
    """
    # yfile = path+'/12-14-17_broadband_src_MZI/interferogram_'+signal+'_v1.txt'
    # yf = pd.read_csv(yfile, sep='\t', usecols=[0,1])
    # yval_train, OPL = yf.values[:,1], yf.values[:,0]
    # yfile = path+'/12-14-17_broadband_src_MZI/interferogram_'+signal+'_v2.txt'
    # yf = pd.read_csv(yfile, sep='\t', usecols=[0,1])
    # yval_validate, OPL = yf.values[:,1], yf.values[:,0]


    """ UNCOMMENT area below only if desired spectrum is KNOWN beforehand
    """
    # xfile = str(signal)+'.CSV'
    # xf = pd.read_csv(xfile, header=30)
    # xval_train, xwl = xf.values[:,1], xf.values[:,0]
    # xwl = np.array([x - 0.7 for x in xwl]) #0.7nm offset between uncalibrated & calibrated OSA
    # wavelengths = np.linspace(1550, 1570, A1.shape[1])
    # x_real = np.interp(wavelengths, xwl, xval_train)
    # y_real = np.dot(A1, x_real)

    """BELOW IS PUTTING IN SIMULATION INFORMATION (no noise at all, only see the accuracy in reconstruction)
    """
    xfile = str(signal) + '.CSV'
    xf = pd.read_csv(xfile, header=30)
    xval_train, xwl = xf.values[:, 1], xf.values[:, 0]
    xwl = np.array([x - 0.7 for x in xwl])  # 0.7nm offset between uncalibrated & calibrated OSA
    wavelengths = np.linspace(1550, 1570, A1.shape[1])
    x_real = np.interp(wavelengths, xwl, xval_train)
    #OPL = path_lengths[sort_order]
    yval_train = np.dot(A1, x_real)
    yval_validate = np.dot(A2, x_real)

    """Adding in noise to the signals yval_train and yval_validate
    """
    noise_level = 0.03
    noise1 = np.random.normal(loc=0.0, scale=np.mean(yval_train)*noise_level, size=A1.shape[0])
    yval_train = np.add(yval_train, noise1)
    noise2 = np.random.normal(loc=0.0, scale=np.mean(yval_train)*noise_level, size=A1.shape[0])
    yval_validate = np.add(yval_validate, noise2)

    """ Pseudo-inverse method (for reference)
    """
    Ainv = np.linalg.pinv(A1)
    x_pinv_train = np.dot(Ainv, yval_train)
    x_pinv_validate = np.dot(Ainv, yval_validate)

    """ Begin ELASTIC-D1 parameter search
    """
    Dsize = len(A1[0])
    D1 = np.zeros((Dsize+1, Dsize))
    for i in range(Dsize):
        D1[i][i] = 1
        D1[i+1][i] = -1

    def elastic_D1_smoothing(l1, l2, l3, A, D1, Dsize, yval):
        """ p = 2(AT.A + l2.I + l3.DT.D) """
        p = 2*(np.dot(np.transpose(A), A) + l2*np.identity(Dsize) + l3*np.dot(np.transpose(D1), D1))
        """ Q = (l1.vec1 - 2 AT.y) """
        Q = l1*np.ones(Dsize) - 2*np.dot(np.transpose(A), yval)
        """ G = -I """
        G = -1*np.identity(Dsize)
        """ h = zero-vector """
        h = np.zeros(Dsize)
        sol = solvers.qp(matrix(p.tolist()), matrix(Q.tolist()),  matrix(G.tolist()), matrix(h.tolist()))
        x = np.transpose(np.array(sol['x']))[0]
        return x

    """ dl, d2, d3 specify the size of the hyperparameter search space
    BE CAREFUL! The total time will be proportional to d1*d2*d3
    Choose 12x12x12 for a search space similar to what is used in the paper
    """
    d1, d2, d3 = 4, 4, 4 # 10, 10, 10 #4,4,4  #12,12,12

    """ Choose a logarithmic range suitably large
    """
    l1_list = np.logspace(-2, 4, d1)
    l2_list = np.logspace(-3, 5, d2)
    l3_list = np.logspace(2, 7, d3)
    r2_list = []
    mv = -1E10
    l1max, l2max, l3max = 0.0, 0.0, 0.0
    time_list = []
    print("Running sweep...")
    for i in range(len(l1_list)):
        l1 = l1_list[i]
        r2_list.append([])
        print("l1="+str(l1))
        print("..progress="+str(100.0*i/len(l1_list))+"%")
        for l2 in l2_list:
            r2_list[-1].append([])
            print(".")
            for l3 in l3_list:
                starttime = timeit.default_timer()
                x = elastic_D1_smoothing(l1, l2, l3, A1, D1, Dsize, yval_train)
                if np.isfinite(x.all()) and max(x) !=0.0:
                    score = r2_score(yval_validate, np.dot(A2, x))
                    r2_list[-1][-1].append(score)
                    if score > mv:
                        mv = score
                        l1max = l1
                        l2max = l2
                        l3max = l3
#                print("compute-time = "+str(timeit.default_timer() - starttime))
                time_list.append(timeit.default_timer() - starttime)
    print("Average time per measurement = "+str(np.average(time_list)))

    print("l1 max = "+str(l1max))
    print("l2 max = "+str(l2max))
    print("l3 max = "+str(l3max))

    xmax = elastic_D1_smoothing(l1max, l2max, l3max, A1, D1, Dsize, yval_train)

    """ Comment below if the spectrum is *not* known beforehand
    """
    # x_real = np.linalg.norm(x_real)*xmax  # CHANGED THIS FROM normalize_vector(x_real, xmax)
    r2_score_ed1 = r2_score(x_real/max(x_real), xmax/max(xmax))
    # l2_score_ed1 = np.linalg.norm((x_real/max(x_real)) - (xmax/max(xmax)), 2)
    print("r^2 result for D1 : %f" % r2_score_ed1)
    # print("L2-norm result for D1: %f" % l2_score_ed1)
    # r2_score_pinv = r2_score(x_real/max(x_real), x_pinv_validate/max(x_pinv_validate))
    # l2_score_pinv = np.linalg.norm((x_real/max(x_real)) - (x_pinv_validate/max(x_pinv_validate)), 2)
    # print("r^2 result for PINV : %f" % r2_score_pinv)
    # print("L2-norm result for PINV: %f" % l2_score_pinv)

    """ Save data to a .csv file if you'd like
    """
#    with open(savedfilename, 'as') as wf:
#        wf.write(str(signal)+', '+str(r2_score_ed1)+', '+str(l2_score_ed1)+', '+str(l1max)+', '+str(l2max)+', '+str(l3max)+'\n')
#        wf.close()


"""============================================================================
Plot overlapping one another
============================================================================"""
plt.figure(figsize=(6,4))

""" Comment first plot below if spectrum is *not* known beforehand
"""
plt.plot(wavelengths, x_real/np.max(x_real), 'k-', label="OSA", linewidth=2.0)
plt.xticks(fontsize=ticksize)
plt.yticks(fontsize=ticksize)
plt.xlabel("Wavelength [nm]", fontsize=labelsize)
plt.ylabel("Intensity [a.u.]", fontsize=labelsize)
maxy, miny = 1.0, 0.0
plt.ylim([miny - 0.05*abs(maxy-miny), maxy + 0.05*abs(maxy-miny)])

plt.plot(wavelengths, xmax/np.max(xmax), 'r-', label="64-ch dFT", linewidth=2.0)
plt.ylabel("Intensity [a.u.]", fontsize=labelsize)
plt.xticks(fontsize=ticksize)
plt.yticks(fontsize=ticksize)
plt.legend(loc='best')
maxy, miny = 1.0, 0.0
plt.ylim([miny - 0.05*abs(maxy-miny), maxy + 0.05*abs(maxy-miny)])

#plt.plot(wavelengths, x_pinv_train/np.max(x_pinv_train), 'b-', label="Pseudoinverse", linewidth=2.0)
#plt.xticks(fontsize=ticksize)
#plt.yticks(fontsize=ticksize)
#plt.xlabel("Wavelength [nm]", fontsize=labelsize)
#plt.ylabel("Intensity [a.u.]", fontsize=labelsize)
#maxy, miny = 1.0, 0.0
#plt.ylim([miny - 0.05*abs(maxy-miny), maxy + 0.05*abs(maxy-miny)])

plt.legend(loc='best', fontsize=12)
plt.tight_layout()


plt.show()
