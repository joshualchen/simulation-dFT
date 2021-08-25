import os
import jcamp
'''jcamp package can be found here: https://pypi.org/project/jcamp/'''
import numpy as np
import math
from matplotlib import pyplot as plt

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def broadSpectrumGauss():
    directory = 'infrared-database'
    y_running_sum = np.zeros((56400,))
    for filename in os.listdir(directory):
        if filename.endswith(".jdx"):
            jcamp_obj = jcamp.JCAMP_reader(directory + "/" + filename)
            x_vals_invcm = jcamp_obj['x']  # units in cm^-1
            y_vals = np.flip(jcamp_obj['y'])[:56400]  # relative units for our purposes
            x_vals = np.flip([2*math.pi/x*10**7 for x in x_vals_invcm])[:56400]
            y_running_sum = np.add(y_running_sum, y_vals)

    intensities = y_running_sum[35000:55001] / 8  # the characterization matrix will be 20001 long. /8 just to scale
    wavelengths = np.linspace(1200, 1700, 20001)
    gauss = np.array([gaussian(x, 1550, 100) for x in wavelengths]) * max(intensities) / 1.3
    intensities_w_gauss = np.add(intensities, gauss)
    return intensities_w_gauss, wavelengths

def broadSpectrum():
    directory = 'infrared-database'
    y_running_sum = np.zeros((56400,))
    for filename in os.listdir(directory):
        if filename.endswith(".jdx"):
            jcamp_obj = jcamp.JCAMP_reader(directory + "/" + filename)
            x_vals_invcm = jcamp_obj['x']  # units in cm^-1
            y_vals = np.flip(jcamp_obj['y'])[:56400]  # relative units for our purposes
            x_vals = np.flip([2*math.pi/x*10**7 for x in x_vals_invcm])[:56400]
            y_running_sum = np.add(y_running_sum, y_vals)

    intensities = y_running_sum[35000:55001] / 8  # the characterization matrix will be 20001 long. /8 just to scale
    wavelengths = np.linspace(1200, 1700, 20001)
    return intensities, wavelengths

