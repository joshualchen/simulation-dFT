#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 14:15:46 2017

@author: joshualchen
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
from create_A_matrix import fast_createA
solvers.options['show_progress']=False

font = {'weight' : 'normal',
        'size'   : 14}
labelsize = 16
ticksize = 14
matplotlib.rc('font', **font)

"""Setting parameters
"""
num_points = 801
min_wavelength = 1550
max_wavelength = 1570
num_trials = 7  # number of times to test each parameter combo (since randomized noise)
pickle_file = "stored/809_broadband_dL_numswitches_sL_vs_R2mean_R2std.pickle"  # where to store info from grid search

""" Initializing
"""
dL_values = []
num_switches_values = []
sL_values = []
R2_mean_values = []
R2_std_values = []
reconstructed_values = []

# uncomment the below if you hope to add to an existing pickle file
# with open(pickle_file, 'rb') as f:
#     [dL_values, R2_values, reconstructed_values, x_real] = pickle.load(f)

"""Load mystery spectrum once, first
"""
signal = "MZI1"
xfile = str(signal) + '.CSV'
xf = pd.read_csv(xfile, header=30)
xval_train, xwl = xf.values[:, 1], xf.values[:, 0]
xwl = np.array([x - 0.7 for x in xwl])  # 0.7nm offset between uncalibrated & calibrated OSA

# start grid search sweep
for dL in np.linspace(5, 55, 6, dtype=int):
    for num_switches in np.linspace(6, 10, 5, dtype=int):
        print("Creating A Matrix...")
        start = time.time()
        A1 = fast_createA(num_points, min_wavelength, max_wavelength, num_switches, dL, 0)
        end = time.time()
        print("Time to create A matrix: " + str(end - start) + " s")
        wavelengths = np.linspace(min_wavelength, max_wavelength, A1.shape[1])
        x_real = np.interp(wavelengths, xwl, xval_train)
        all_R2_vals = []
        for i in range(num_trials):
            start = time.time()

            yval_train = np.dot(A1, x_real)

            """Adding in noise to the signals yval_train and yval_validate
            """
            noise_level = 0.03
            noise1 = np.random.normal(loc=0.0, scale=np.mean(yval_train)*noise_level, size=A1.shape[0])
            yval_train = np.add(yval_train, noise1)

            """ Pseudo-inverse method (for reference)
            """
            Ainv = np.linalg.pinv(A1)
            x_pinv_train = np.dot(Ainv, yval_train)

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
                            score = r2_score(x, x_real)
                            r2_list[-1][-1].append(score)
                            if score > mv:
                                mv = score
                                l1max = l1
                                l2max = l2
                                l3max = l3
        #                print("compute-time = "+str(timeit.default_timer() - starttime))
                        time_list.append(timeit.default_timer() - starttime)

            xmax = elastic_D1_smoothing(l1max, l2max, l3max, A1, D1, Dsize, yval_train)

            """ Comment below if the spectrum is *not* known beforehand
            """
            r2_score_ed1 = r2_score(x_real/max(x_real), xmax/max(xmax))
            print("r^2 result for D1 : %f" % r2_score_ed1)
            all_R2_vals.append(r2_score_ed1)

        R2_mean_values.append(np.mean(all_R2_vals))
        R2_std_values.append(np.std(all_R2_vals))
        num_switches_values.append(num_switches)
        dL_values.append(dL)
        sL_values.append(0)

# store data in pickle file
os.makedirs(os.path.dirname(pickle_file), exist_ok=True)
with open(pickle_file, 'wb') as f:
    pickle.dump([dL_values, num_switches_values, sL_values, R2_mean_values, R2_std_values], f)
