from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import scipy.interpolate
import pickle

pickle_file = "stored/811_broadband_dL_numswitches_sL_vs_R2mean_R2std.pickle"

with open(pickle_file, 'rb') as f:
    [dL_values, num_switches_values, sL_values, R2_mean_values, R2_std_values] = pickle.load(f)

# Here, dL and num_switches are linear, sL is exponential
# x axis are the num_switches values, y axis are the dL values

sL = 0
indices = np.where(np.array(sL_values) == sL)[0]
#indices = np.linspace(1, 93, 24, dtype=int)

# Define the values to be plotted
num_switches_values_sL = np.array(num_switches_values)[indices]
dL_values_sL = np.array(dL_values)[indices]
sL_values_sL = np.array(sL_values)[indices]
R2_mean_values_sL = np.array(R2_mean_values)[indices]
R2_std_values_sL = np.array(R2_std_values)[indices]

xi, yi = np.linspace(num_switches_values_sL.min(), num_switches_values_sL.max(), 100), np.linspace(dL_values_sL.min(), dL_values_sL.max(), 100)
xi, yi = np.meshgrid(xi, yi)

rbf_mean = scipy.interpolate.Rbf(num_switches_values_sL, dL_values_sL, R2_mean_values_sL, function='linear')
rbf_std = scipy.interpolate.Rbf(num_switches_values_sL, dL_values_sL, R2_std_values_sL, function='linear')
zi_mean = rbf_mean(xi, yi)
zi_std = rbf_std(xi, yi)

# norm below may potentially be good as well
#norm = mcolors.TwoSlopeNorm(vmin=R2_std_values_sL.min(), vcenter=0.01, vmax=R2_std_values_sL.max())

plt.figure("Mean", figsize=[5, 5])
norm = mcolors.PowerNorm(gamma=0.5)
cmap = 'nipy_spectral_r'
#plt.imshow(zi_mean.max()-zi_mean, norm=norm,
#           origin='lower', extent=[num_switches_values_sL.min(), num_switches_values_sL.max(), dL_values_sL.min(), dL_values_sL.max()], aspect='auto', cmap=cmap)
plt.scatter(num_switches_values_sL, dL_values_sL, c=zi_mean.max()-R2_mean_values_sL, cmap=cmap, norm=norm)
for X, Y, Z in zip(num_switches_values_sL, dL_values_sL, R2_mean_values_sL):
    plt.annotate('{}'.format(round(Z, 4)), xy=(X,Y), xytext=(15, 5), ha='right', textcoords='offset points')
cb = plt.colorbar()
plt.title("Mean R2 value, plotting num switches vs. dL value")
plt.xlabel("Number of Switches")
plt.ylabel("dL Value (um)")
plt.xlim([3.5, 12.5])
plt.ylim([0, 35])

plt.figure("Std", figsize=[5, 5])
cmap = 'nipy_spectral'
#plt.imshow(zi_std, origin='lower', #norm=norm,
#           extent=[num_switches_values_sL.min(), num_switches_values_sL.max(), dL_values_sL.min(), dL_values_sL.max()], aspect='auto', cmap=cmap)
plt.scatter(num_switches_values_sL, dL_values_sL, c=R2_std_values_sL, cmap=cmap)#, norm=norm)
for X, Y, Z in zip(num_switches_values_sL, dL_values_sL, R2_std_values_sL):
    plt.annotate('{}'.format(round(Z, 4)), xy=(X,Y), xytext=(15, 5), ha='right', textcoords='offset points')
cb = plt.colorbar()
plt.title("Std of R2 value, plotting num switches vs. dL value")
plt.xlim([3.5, 12.5])
plt.ylim([0, 35])

plt.show()
