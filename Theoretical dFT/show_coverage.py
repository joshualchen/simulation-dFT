from matplotlib import pyplot as plt
import numpy as np

# This python file takes in parameters that define the characterization matrix, and returns a plot showing the range of
# frequencies that is covered, and where it lies in terms of the bounds we have


def calculateVals(num_switches, dL, min_wavelength, max_wavelength, num_points, n=1.82):
    mid_wavelength = (min_wavelength + max_wavelength) / 2  # in nm
    bottom_OPD = mid_wavelength ** 2 / (2 * n * (max_wavelength - min_wavelength)) / 1000  # in um
    top_OPD = (mid_wavelength ** 2 * num_points) / (4*1000*n*(max_wavelength-min_wavelength))  # in um
    OPD_vals = np.array(range(2 ** num_switches)) * dL  # in um
    return bottom_OPD, top_OPD, OPD_vals


if __name__ == "__main__":
    num_switches = 6
    dL = 26
    min_wavelength = 1550
    max_wavelength = 1570
    num_points = 801

    bottom_OPD, top_OPD, OPD_vals = calculateVals(num_switches, dL, min_wavelength, max_wavelength, num_points)
    print("bottom: " + str(bottom_OPD))
    print("top: " + str(top_OPD))

    y = np.ones(np.shape(OPD_vals))
    plt.plot(OPD_vals, y, '|', ms=40, color='b', label='OPD values present', linewidth='0.5')
    plt.plot(bottom_OPD, 1., '|', ms=40, color='r', label='OPD that gives FSR of 2x the bandwidth', linewidth='0.5')
    plt.plot(top_OPD, 1., '|', ms=40, color='g', label='OPD that gives FSR of 4x the resolution', linewidth='0.5')
    plt.xlabel("OPD Value (um)")
    plt.legend()
    plt.show()
