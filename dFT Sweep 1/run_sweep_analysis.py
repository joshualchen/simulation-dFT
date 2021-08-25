import numpy as np
from matplotlib import pyplot as plt
import pickle

pickle_file = 'stored/711_64_[fsrmean_fsrmiddle_T]_pathlengths_statetable_stateheader.pickle'
# MUST CHANGE PICKLE LOAD ORDER IN LINE 9 if variables stored in different order

with open(pickle_file, 'rb') as f:
    [result_fsr_mean, result_fsr_middle, result_T], path_lengths, state_table, state_header = pickle.load(f)

FSR_middle = result_fsr_middle['FSR_middle']  # FSR_middle is FSR taken at ~1550 nm (the middle point)
FSR_mean = result_fsr_mean['FSR_mean']  # FSR_mean is the mean of all FSR's
together = np.concatenate(([path_lengths], [FSR_middle]), axis=0)  # pairing the FSR_middle values with respective paths
sorted = together[:, together[0].argsort()]  # sorting the path lengths from shortest to largest with FSR_middle
transmission_spectra = np.transpose(result_T['mode 1 transmission'])  # stores all of the spectra from each run

# plot FSR vs path lengths
plt.clf()
plt.figure(1)
plt.scatter(sorted[0], sorted[1])
plt.xlabel("Optical Path Length Difference (um)")
plt.ylabel("FSR @ ~1550 nm (m)")
plt.title("FSR vs. Path Length for 64 states")

# plot the transmission spectrum of a specific plot
plt.figure(2)
index = 1 - 1
axis_min = 1509.41
axis_max = 1589.43
transmission = np.square(np.absolute(transmission_spectra[index]))
wavelengths = np.linspace(axis_min, axis_max, num=len(transmission))
plt.scatter(wavelengths, transmission)
plt.ylabel("Abs^2 Intensity")
plt.xlabel("Wavelength (nm)")
plt.show()
