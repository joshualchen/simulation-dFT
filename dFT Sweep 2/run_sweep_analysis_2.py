import numpy as np
from matplotlib import pyplot as plt
import pickle
from matplotlib.widgets import TextBox
import sys

pickle_file = 'stored/802_16_[all_results]_pathlengths_statetable_stateheader.pickle'
# MUST CHANGE PICKLE LOAD ORDER IN LINE 11 if variables stored in different order

with open(pickle_file, 'rb') as f:
    [all_results, path_lengths, state_table, state_header] = pickle.load(f)

# Structure of all_results: (each line going down represents a deeper structure)
# At the top: "all_results" (*dict* with the name of the sweep as the keys)
# Next: a single "sweep" (*dict* with the name of each result in the sweep as the keys)
# Next: the specific "result" (*dict* with all of the Lumerical outputs. One has the desired info for each result)

# Below, we use this current "all_results" structure and extract the info for each result for a few example plots
# If there are additional results, edit this script to plot what's desired

# Plotting all of the FSR vs. path length plots
sweep_names = list(all_results.keys())
# Plotting a figure for each FSR vs. path length plot
for i in range(len(sweep_names)):
    sweep_dict = all_results[sweep_names[i]]  # getting each sweep dictionary
    result_FSR_middle = sweep_dict['fsr_middle']
    FSR_middle = result_FSR_middle['FSR_middle']
    result_FSR_mean = sweep_dict['fsr_mean']
    FSR_mean = result_FSR_mean['FSR_mean']
    together = np.concatenate(([path_lengths], [FSR_middle]), axis=0)  # placing them together to be sorted
    sorted = together[:, together[0].argsort()]  # sorting the points all together
    plt.figure(sweep_names[i] + "FSR vs. path length")
    plt.scatter(sorted[0], sorted[1])
    plt.xlabel("Optical Path Length Difference (um)")
    plt.ylabel("FSR @ ~1550 nm (m)")
    plt.title("FSR vs. Path Length for Sweep "+sweep_names[i])

#plt.show()

# plotting the transmission of a select sweep, and then the path length, using plotly
# first, create the data structure of a dict of sweeps, within is a 2-column list of path length with transmission
transmission_dict = {}
sweep_names = list(all_results.keys())
for i in range(len(sweep_names)):
    # path_lengths is a 1x16 numpy array
    sweep_dict = all_results[sweep_names[i]]
    T_dict = sweep_dict['T']
    wavelengths_flipped = T_dict['wavelength']
    wavelengths = np.transpose(wavelengths_flipped)[0] * 10 ** 9
    transmission_spectra_flipped = T_dict['mode 1 transmission']
    transmission_spectra = np.transpose(transmission_spectra_flipped)  # now, put an index in to get the transmission
    combined_list = []
    for j in range(len(path_lengths)):
        pair = [path_lengths[j], transmission_spectra[j]]
        combined_list.append(pair)
    combined_np = np.array(combined_list, dtype=object)
    sorted = combined_np[combined_np[:, 0].argsort()]  # first value in array is path length, second is np transmission
    sweep_info = dict()
    sweep_info['wavelength'] = wavelengths
    sweep_info['paths_transmission'] = sorted
    transmission_dict[sweep_names[i]] = sweep_info  # for each sweep key, now is stored a 16x2 shape

# ALL PLOTS IN TRANSMISSION_DICT
# now, plotting with textboxes in matplotlib:
fig, ax = plt.subplots()
fig.canvas.set_window_title('View Transmissions')
fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)
fig.subplots_adjust(bottom=0.25)
new_path_lengths = np.copy(path_lengths)
new_path_lengths.sort()
text_1 = str(list(transmission_dict.keys())[0])
text_2 = str(new_path_lengths[0])
find_index = list(transmission_dict[text_1]['paths_transmission'][:, 0]).index(int(text_2))
y = np.square(np.absolute(transmission_dict[text_1]['paths_transmission'][find_index][1]))
x = transmission_dict[text_1]['wavelength']
l, = ax.plot(x, y)
ax.set_ylabel("Abs^2 Intensity")
ax.set_xlabel("Wavelength (nm)")
ax.title.set_text('Transmission for "'+text_1+'", path difference of '+text_2)
print("Sweep names: " + str(sweep_names).strip("[]"))
print("Path differences: " + str(list(new_path_lengths)).strip("[]"))


def on_press(event):
    sys.stdout.flush()
    global text_1
    global text_2
    index_1 = sweep_names.index(text_1)
    index_2 = list(new_path_lengths).index(int(text_2))
    if event.key == 'left' or event.key == 'right' or event.key == 'up' or event.key == 'down':
        if event.key == 'left':
            if index_2 <= 0:
                text_2 = str(new_path_lengths[0])
            else:
                text_2 = str(new_path_lengths[index_2 - 1])
        elif event.key == 'right':
            if index_2 >= len(new_path_lengths) - 1:
                text_2 = str(new_path_lengths[len(new_path_lengths) - 1])
            else:
                text_2 = str(new_path_lengths[index_2 + 1])
        elif event.key == 'up':
            if index_1 >= len(sweep_names) - 1:
                text_1 = sweep_names[len(sweep_names) - 1]
            else:
                text_1 = sweep_names[index_1 + 1]
        elif event.key == 'down':
            if index_1 <= 0:
                text_1 = sweep_names[0]
            else:
                text_1 = sweep_names[index_1 - 1]
        find_index = list(transmission_dict[text_1]['paths_transmission'][:, 0]).index(int(text_2))
        y = np.square(np.absolute(transmission_dict[text_1]['paths_transmission'][find_index][1]))
        l.set_ydata(y)
        text_box_1.set_val(text_1)
        text_box_2.set_val(text_2)
        ax.relim()
        ax.autoscale_view()
        plt.draw()
        ax.title.set_text('Transmission for "' + text_1 + '", path difference of ' + text_2)
        print("Plotted: " + text_1 + ", " + text_2)


def submit(text1, text2):
    global text_1
    global text_2
    text_1 = text1
    text_2 = text2
    text_box_1.text_disp.set_color('black')
    text_box_2.text_disp.set_color('black')
    ax.title.set_color('black')
    try:
        find_index = list(transmission_dict[text_1]['paths_transmission'][:, 0]).index(int(text_2))
        y = np.square(np.absolute(transmission_dict[text_1]['paths_transmission'][find_index][1]))
        l.set_ydata(y)
        ax.relim()
        ax.autoscale_view()
        plt.draw()
        ax.title.set_text('Transmission for "' + text_1 + '", path difference of ' + text_2)
        print("Plotted: " + text_1 + ", " + text_2)
    except ValueError:
        ax.title.set_text("Error: Path difference not found.")
        ax.title.set_color('red')
        text_box_2.text_disp.set_color('red')
        print("Error: Path difference not found.")
        print("Path differences: " + str(list(new_path_lengths)).strip("[]"))
    except KeyError:
        ax.title.set_text("Error: Path difference not found.")
        ax.title.set_color('red')
        text_box_1.text_disp.set_color('red')
        print("Error: Sweep name not found.")
        print("Sweep names: " + str(sweep_names).strip("[]"))


axbox_1 = fig.add_axes([0.2, 0.05, 0.3, 0.075])
axbox_2 = fig.add_axes([0.65, 0.05, 0.25, 0.075])
text_box_1 = TextBox(axbox_1, "Sweep Name: ", initial=text_1)
text_box_1.on_submit(lambda text: submit(text, text_2))
text_box_2 = TextBox(axbox_2, "Path Diff: ", initial=text_2)
text_box_2.on_submit(lambda text: submit(text_1, text))
fig.canvas.mpl_connect('key_press_event', on_press)
plt.show()
