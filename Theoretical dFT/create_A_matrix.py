import numpy as np
import math
import time

# create A matrix just from number of transmissions, num points wanted, dL value

num_switches = 6
min_x = 1550
max_x = 1570
num_points = 801
dL = 40
sL = 0

# x_arr: the array of x values to sample, in nanometers
# a_arr: a is the multiple of dL for this line (OPD = a * dL)
# dL: the unit length difference, in microns
# sL: shift, where L = a * dL + sL
# n: refractive index
# s_arr: the array of phase shift values s that is added to one arm
# D: distance at measurement, in microns


def createA(num_points, min_x, max_x, num_switches, dL, sL):
    def createMatrix(x_arr, a_arr, dL, sL, n, s_arr, D):
        def oneLine(x_arr, a, dL, sL, n, s, D):
            E1r = lambda x: math.cos(2 * math.pi * D / (x / 1000))
            E1i = lambda x: math.sin(2 * math.pi * D / (x / 1000))
            E2r = lambda x: math.cos(2 * math.pi * (D + n * a * dL + sL) / (x / 1000) + s)
            E2i = lambda x: math.sin(2 * math.pi * (D + n * a * dL + sL) / (x / 1000) + s)
            I = lambda x: (E1r(x) + E2r(x)) ** 2 + (E1i(x) + E2i(x)) ** 2
            return np.array([I(x) for x in x_arr])
        return np.concatenate(np.array([np.array([oneLine(x_arr, a, dL, sL, n, s, D) for a in a_arr]) for s in s_arr]), axis=0)
    n = 1.82
    D = 1
    s_arr = [0, math.pi/2]
    # dL = 33.43  # normally, "optimal(?)" dL may be calculated via bounds, but we set for now
    x_arr = np.linspace(min_x, max_x, num_points)
    a_arr = range(2 ** num_switches)
    return createMatrix(x_arr, a_arr, dL, sL, n, s_arr, D)


# equations are simplified rather than written out fully
def fast_createA(num_points, min_x, max_x, num_switches, dL, sL):
    def createMatrix(x_arr, a_arr, dL, sL, n, s_arr):
        def oneLine(x_arr, a, dL, sL, n, s):
            I = lambda x: 2 * math.cos(2 * math.pi * (n * a * dL + sL) / (x / 1000) + s) + 2
            return np.array([I(x) for x in x_arr])
        return np.concatenate(np.array([np.array([oneLine(x_arr, a, dL, sL, n, s) for a in a_arr]) for s in s_arr]), axis=0)
    n = 1.82
    s_arr = [0, math.pi/2]
    # dL = 33.43  # normally, "optimal(?)" dL may be calculated via bounds, but we set for now
    x_arr = np.linspace(min_x, max_x, num_points)
    a_arr = range(2 ** num_switches)
    return createMatrix(x_arr, a_arr, dL, sL, n, s_arr)


if __name__ == "__main__":
    print("Creating A Matrix...")
    start = time.time()
    A1 = fast_createA(num_points, min_x, max_x, num_switches, dL, sL)
    end = time.time()
    print("Time to create A matrix: " + str(end - start) + " s")
    print("A1 is "+str(A1.shape[0])+"x"+str(A1.shape[1])+". Note: A1 contains the sine counterpart (additional phase shifter on one arm of MZI).")
