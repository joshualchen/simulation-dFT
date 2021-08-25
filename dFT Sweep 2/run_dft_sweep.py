import csv
import numpy as np
import pickle
import time
import os
import shutil
import sys
sys.path.append("C:/Program Files/Lumerical/v212/api/python")
import lumapi

# VERSION 2: USING FREQ-DEPENDENT PHASE SHIFTERS INSTEAD OF VOLTAGE
# Main changes:
#   - switch_table_file holds strings (e.g. "pi/2", "3pi/2") instead of number voltage values
#   - must map input strings to text files that indicate frequency vs. phase shift
#   - according to algorithm parameters, all frequency vs. phase shift txt files will be created before the sweep
#   - sweep will include a *new layer* that sweeps "center wavelengths", due to frequency dependent phase shifts
#   - calculation for numPoints() changes to address the new layer of center wavelengths

# Switch files below to save the output results from sweep
simulation_file = "equinor_3stage.icp"  # icp file referenced by switch_table_file
ONA_name = "ONA_1"  # what is the ONA object in the simulation file called?
switch_table_file = "switch_table_64.csv"  # csv file formatted according to guidelines
pickle_file = 'stored/802_64_[all_results]_pathlengths_statetable_stateheader.pickle'  # stored variables location
# Calculate path differences using top_path & bottom_path & switch_paths in function calculatePathLengths()
top_path = [0, 1, 2]  # these array elements refer to indices of switch_paths & other lists defined in calculatePathLengths()
bottom_path = [5, 6, 7]  # DEPENDING ON CIRCUIT DESIGN, PATH DIFFERENCE MAY BE CALCULATED DIFFERENTLY
switch_paths = np.array([[40, 0], [160, 0], [640, 0], [2560, 0], [10240, 0], [40, 120], [160, 480], [640, 1920],
                         [2560, 7680], [10240, 30720]])  # defines path lengths for each switch: [path up, path down]
# Indicate parameters for wavelength segments given the frequency dependent phase shifts
min_wavelength = 1675  # min and max wavelengths indicate the edges of the broadband
max_wavelength = 1800
num_segments = 2  # the number of segments to split the broadband into
phase_shifts = [["pi2", "pi", "3pi2"], [1.570796, 3.1415927, 4.712389]]  # name of the phase shift, then the value
# Choose the results in results_dict to be collected by the sweep
results_dict = [
    # before the sweep can be run, fsr_mean & fsr_middle must be extracted by code from the root analysis script
    # above is done automatically at the start of runSweep

    # LUMERICAL RESULT DICT KEYS FORMAT:
    # Name - chosen name for Python reference
    # Result - must refer to a property/value within Lumerical

    {
        "Name": "fsr_mean",
        "Result": "::Root Element::FSR_mean"
    },

    {
        "Name": "fsr_middle",
        "Result": "::Root Element::FSR_middle"
    },
    # make sure that references to the components are correct (e.g. ONA_3 vs. ONA_1)
    {
        "Name": "T",
        "Result": "::Root Element::"+ONA_name+"::input 1/mode 1/transmission",
        "Estimation": True
     }
]

# QUICK DEFINITIONS/EXPLANATIONS:
# switch_table is the initial user input table that defines the switches and the states for each switch
# state_table is created from switch_table and lists every possible state for the parameter sweep, each uniquely associated with a dL
# state_header indicates the columns of state_table

# input switch_table_file (csv) and convert to state_table
def createStateTable(switch_table_file):
    # state_header: provides the switch and property strings associated with the columns of state_table
    # state_table: has (total # of switch combinations) # of rows, (total # property switch combos) # of columns

    # import using csv reader
    switch_table = list(csv.reader(open(switch_table_file)))

    # listing all of the switches to put together the columns of state_table
    switch_list = [i[0] for i in switch_table[1:]]

    # create and fill state_header 2d string array (2 cols, total # property switch combos rows)
    state_header = []
    for i in range(len(switch_list)):
        # must use string manipulation to separate csv inputs into list elements
        property_vals = [s.strip() for s in switch_table[i+1][1].strip().strip("[]").split(",")]
        # listing all switch name and switch property combinations together in columns of two rows
        for m in range(len(property_vals)):
            state_header.append([switch_list[i], property_vals[m]])

    # fillTable builds state_table from switch_table using recursion
    def fillTable(existing_matrix, switches_left, switch_table):
        if not switches_left:
            # once all switches have been incorporated, return the accumulated state_table
            return existing_matrix
        else:
            # get height of existing structure to see number of times to repeat new column
            height_existing_matrix = existing_matrix.shape[0]
            # choose the next switch to format into the table
            next_switch = switches_left[0]
            # find the index of the switching block in the switch_table
            helper_index = [i[0] for i in switch_table[1:]].index(next_switch)
            # list the states for this switching block
            states = [a for a in switch_table[helper_index+1][2:] if a.strip()]
            # the number of states for this switching block
            num_states = len(states)
            # convert states to np array format
            np_states = np.array([[str(s.strip().strip('""')) for s in state.strip().strip("[]").split(",")] for state in states])
            # create the form of the values to be appended from the states as np array in new matrix
            new_states = np.concatenate([np_states]*height_existing_matrix, axis=0)
            # must take existing states and repeat each row separately
            repeated_existing = np.repeat(existing_matrix, repeats=num_states, axis=0)
            # combine the new states with existing states to form the new matrix
            new_matrix = np.concatenate((repeated_existing, new_states), axis=1)
            # recurse to continue down properties_left and build the new matrix
            return fillTable(new_matrix, switches_left[1:], switch_table)

    # call fillTable to return state_table, existing matrix is empty initially
    state_table = fillTable(np.empty((1, 0)), switch_list, switch_table)

    return state_table, state_header


# calculate the number of points needed for the sweep: smallest wavelength to be measured, largest path difference
def numPoints(min_wavelength, max_wavelength, group_index, top_path, bottom_path, switch_paths):
    # min_wavelength/max_wavelenght in nm, switch_paths in um, theoretical FSR in nm
    smallest_top_path = 0
    for switch in top_path:
        smallest_top_path += min(switch_paths[switch])
    largest_bottom_path = 0
    for switch in bottom_path:
        largest_bottom_path += max(switch_paths[switch])
    # calculating largest path difference; correlates to smallest FSR
    max_path_diff = largest_bottom_path - smallest_top_path  # paths in um
    # using the min wavelength because FSR to become a minimum
    theoretical_FSR = min_wavelength ** 2 / (group_index * max_path_diff) * 10 ** -3  # FSR in nm
    # number of points necessary per FSR
    points_to_FSR = 25
    # the necessary number of points to span the entire band, at least 10 points per FSR
    num_points = round((max_wavelength - min_wavelength) / theoretical_FSR * points_to_FSR)
    print("Number of points in ONA: "+str(num_points))
    return num_points


# create the frequency vs. phase shift text files for each center wavelength
def createPhaseFiles(min_wavelength, max_wavelength, num_segments, phase_shifts):
    # delete phase_shift_files folder, and then repopulate
    try:  # try to remove, if no directory exists then print in terminal
        shutil.rmtree("phase_shift_files")
    except OSError as e:
        print("Delete phase_shift_files: no directory exists")
    # define start frequency and stop frequency to create frequency array
    start_freq = 2.998 * 10 ** 17 / max_wavelength
    stop_freq = 2.998 * 10 ** 17 / min_wavelength
    # first, create the center frequencies array
    start_points = np.linspace(min_wavelength, max_wavelength, num=num_segments, endpoint=False)
    mid_points = start_points + (max_wavelength - min_wavelength) / (2 * num_segments)
    inc = 500000000000 # *******CAN WE MAKE THIS VALUE EVEN BIGGER/INTERPOLATE EVEN MORE?? WOULD MAKE READING FILES SO MUCH FASTER
    all_freq_values = np.arange(start_freq, stop_freq + inc, inc)
    # for each center value, we want to create a folder and have the target phase folders within
    for i in range(0, len(mid_points)):
        # define the frequency in nm
        wavelength = round(mid_points[i])
        # for each phase, we want to create a file
        for j in range(0, len(phase_shifts[0])):
            # create the 2d array of freq vs. phase shift
            phase_array = all_freq_values * phase_shifts[1][j] / (2.998 * 10 ** 17 / mid_points[i])
            together = np.concatenate(([all_freq_values], [phase_array]), axis=0)
            together_vert = np.transpose(together)
            # defining filename, including directories
            filename = "phase_shift_files/freq_vs_phase_shift_"+str(wavelength)+"nm/"+str(wavelength)+"_"+phase_shifts[0][j]+"_shift.txt"
            # making sure the files are created
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            # writing to the text file
            with open(filename, 'wb') as f:
                np.savetxt(f, together_vert, delimiter=' ', newline='\n', header='', footer='', comments='# ')
    # return mid_points in units wavelengths so that we can loop through
    return mid_points

# opens Lumerical and executes the sweep, returns the results of the sweep
# inputs: .icp project file, the name of the sweep, the csv file with switch specifications, and the results to measure
def runSweep(project_file, ONA_name, sweep_name, state_table, state_header, results_dict, top_path, bottom_path, switch_paths, min_wavelength, max_wavelength, num_segments, phase_shifts):
    # create interconnect component through lumerical api
    ic = lumapi.INTERCONNECT(filename=project_file, hide=False)
    ic.switchtodesign()

    # adding FSR_mean & FSR_middle to the root element to sweep, must create new property
    ic.select("::Root Element")
    # string set to property "analysis script" is lumerical code to create the FSR_mean & FSR_middle properties
    ic.set("analysis script", 'FSR_data = getresult("::Root Element::'+ONA_name+'","input 1/mode 1/peak/free spectral '
                              + 'range"); ?FSR_mean = mean(FSR_data.getattribute(getattribute(FSR_data))); '
                              + 'setresult("FSR_mean", FSR_mean); '
                              + 'FSR_list = FSR_data.getattribute(getattribute(FSR_data)); '
                              + '?FSR_middle = FSR_list(length(FSR_list)/2); setresult("FSR_middle", FSR_middle);')
    # must set ONA values to min and max wavelength, and then set number of points
    ic.setnamed(ONA_name, "input parameter", "start and stop")
    start_freq = 2.998 * 10 ** 17 / max_wavelength
    stop_freq = 2.998 * 10 ** 17 / min_wavelength
    ic.setnamed(ONA_name, "start frequency", start_freq)
    ic.setnamed(ONA_name, "stop frequency", stop_freq)
    group_index = 2  # May be subject to change, may need to move up for easy access as well
    # calculate the minimum number of points to have good enough FSR resolution
    num_points = numPoints(min_wavelength, max_wavelength, group_index, top_path, bottom_path, switch_paths)
    ic.setnamed(ONA_name, "number of points", 20000)

    # measure the amount of time for one run to predict total run time NEED TO FIX THIS, INACCURATE NOW
    start = time.time()
    ic.run()
    end = time.time()
    run_time_one = end - start

    # create all of the file folders, and pass in mid_points for looping
    mid_points = createPhaseFiles(min_wavelength, max_wavelength, num_segments, phase_shifts)

    # total run time estimate is the number of sweeps to run * number of values in the sweep * time for one run * 2 (ratio)
    run_time_total = 2 * len(mid_points) * state_table.shape[0] * run_time_one
    m, s = divmod(run_time_total, 60)
    h, m = divmod(m, 60)
    print("Running sweep. Crude time estimate: " + str(h) + " hr, " + str(m) + " min")

    # initiate array to store all of the different sweep results
    all_results = {}
    # creating a for loop to create N number of sweeps, one for each middle wavelength
    for i in range(0, len(mid_points)):
        # name of new sweep with center wavelength appended
        segment_sweep_name = sweep_name + "_" + str(round(mid_points[i]))

        # setup sweep
        # delete previous instances of sweep
        ic.deletesweep(segment_sweep_name)
        ic.addsweep(0)
        ic.setsweep("sweep", "name", segment_sweep_name)
        ic.setsweep(segment_sweep_name, "type", "Values")
        ic.setsweep(segment_sweep_name, "number of points", state_table.shape[0])
        print("Sweep created: " + ic.getsweep(segment_sweep_name, "name"))

        # setup sweep parameters according to state_table and state_header
        # add all sweep parameters as dictionaries
        for j in range(len(state_header)):
            object_name = state_header[j][0]
            property_name = state_header[j][1]
            para = {
                "Name": object_name + ", " + property_name,
                "Parameter": "::Root Element::" + object_name + "::" + property_name,
                "Type": "Text"
            }
            # fill the parameter settings for each object/property combo using data from state_table
            for k in range(state_table.shape[0]):
                # define the desired phase_shift from the table
                phase_shift = np.transpose(state_table)[j][k]
                # extract the center_wavelength for this component
                center_wavelength = round(mid_points[i])
                # defined filename according to how it was created
                filename = "'"+os.getcwd()+"/phase_shift_files/freq_vs_phase_shift_"+str(center_wavelength)+"nm/"+str(center_wavelength)+"_"+phase_shift+"_shift.txt'"
                # assign the value to the property
                para["Value_" + str(k + 1)] = filename
            # add the created dictionary as another sweep parameter
            ic.addsweepparameter(segment_sweep_name, para)

        # set up results parameters
        for result in results_dict:
            ic.addsweepresult(segment_sweep_name, result)

        # run the sweep
        ic.runsweep(segment_sweep_name)
        # get sweep results and return in a dict, to be placed in a dict with other sweeps
        results = {}
        for result in results_dict:
            results[result["Name"]] = ic.getsweepresult(segment_sweep_name, result["Name"])

        all_results[segment_sweep_name] = results

    # Structure of all_results: (each line going down represents a deeper structure)
    # At the top: "all_results" (*dict* with the name of the sweep as the keys)
    # Next: a single "sweep" (*dict* with the name of each result in the sweep as the keys)
    # Next: the specific "result" (*dict* with all of the Lumerical outputs. One has the desired info for each result)
    return all_results


# given state_table, calculate path lengths for each switch combo (MUST CHANGE BASED ON CHIP DESIGN)
# currently, the function below calculates path differences for two-path MZI designs with specific voltage values
def calculatePathLengths(state_table, top_path, bottom_path, switch_paths):
    # top_path and bottom_path define the top & bottom switches, with indices corresponding to switch_paths below
    # DEFINED AT TOP OF FILE FOR EASY ACCESS, BUT COMMENTED BELOW:
    # top_path = [0, 1, 2, 3, 4]
    # bottom_path = [5, 6, 7, 8, 9]
    # switch_paths AGAIN, COMMENTED BELOW, BROUGHT TO TOP FOR EASY ACCESS:
    # switch_paths defines path lengths for each switch; first value is path length for up, second is for down
    # switch_paths = np.array([[40, 0], [160, 0], [640, 0], [2560, 0], [10240, 0], [40, 120], [160, 480], [640, 1920],
    #                        [2560, 7680], [10240, 30720]])

    # switch_type defines the type of switch (with indices mirroring switch_paths).
    # 0 means the switch setting defines up or down, 1 means the switch setting defines switch up/down or stay up/down
    switch_type = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    # indicator also mirrors indices of switch_paths & switch_type
    # indicator lists for each switch the setting for "go up" for switch_type of 0, or "stay" for switch_type of 1; i.e. "State 1"
    indicator = [['3pi2','pi2'],['3pi2','pi2'],['3pi2','pi2'],['3pi2','pi2'],['3pi2','pi2'],['3pi2','pi2'],['3pi2','pi2'],['3pi2','pi2'],['3pi2','pi2'],['3pi2','pi2']] # up or stay indicator
    # indicator_not lists the setting for "go down" for switch_type of 0, or "switch" for switch_type of 1; i.e. "State 2"
    indicator_not = [['pi2','pi2'],['pi2','pi2'],['pi2','pi2'],['pi2','pi2'],['pi2','pi2'],['pi2','pi2'],['pi2','pi2'],['pi2','pi2'],['pi2','pi2'],['pi2','pi2']] # down or switch indicator

    # create path length difference array, append a value for each configuration listed in state_table
    path_differences = []
    # for each configuration listed in state_table
    for i in range(state_table.shape[0]):
        # config lists the switch voltage settings in a concatenated list (no distinction between switches)
        config = state_table[i]
        # initializing top length as 0 to build the total top path length by iterating through the top switches
        top_length = 0
        # running_down tracks the current direction the light is traveling (up/down)
        running_down = True  # True means going down, False means going up
        # flattening into a list that we work through each iteration
        running_state_top = [config[j] for j in range(0, str([indicator[i] for i in top_path]).count(",") + 1)]
        starting_len_top = len(running_state_top)
        for j in top_path:
            num_values = len(indicator[j])
            if running_state_top[0:num_values] == indicator[j]:  # it means we're up or staying
                if switch_type[j] == 0:  # this means just go in the direction
                    running_down = False
                    top_length += switch_paths[j][int(running_down)]
                elif switch_type[j] == 1:  # if switch_type is 1
                    # running_down stays the same
                    top_length += switch_paths[j][int(running_down)]
            elif running_state_top[0:num_values] == indicator_not[j]: # if it doesn't match the indicator, as in it's down or switching
                if switch_type[j] == 0:  # this means it just goes down
                    running_down = True
                    top_length += switch_paths[j][int(running_down)]
                elif switch_type[j] == 1:  # this means it's switching
                    running_down = not running_down
                    top_length += switch_paths[j][int(running_down)]
            else:
                # If this error is printed during runtime, then something defined in this function is wrong
                print("ERROR: Check indicator array with state_table array; state mismatch.")
            running_state_top = running_state_top[num_values:]
        # at end of for loop, the top length for one config is held in top_length
        # all we need to do is do the same for the bottom length, and then hold it in bottom_length
        # then subtract top_length from bottom_length, and then append the value into path_lengths

        # Check to see if all values have been accounted for from state_table
        if len(running_state_top) != 0:
            # If this error is printed during runtime, something defined in this function is wrong
            print("ERROR: Not all states have been accounted for in top_path; state mismatch")

        # same process as calculating top_length above
        bottom_length = 0
        running_down = True
        # must accumulate running_state for bottom path from state_table, but only after the values from top path
        running_state_bottom = [config[j] for j in range(starting_len_top, starting_len_top + str([indicator[i] for i in bottom_path]).count(",") + 1)]
        for j in bottom_path:
            num_values = len(indicator[j])
            if running_state_bottom[0:num_values] == indicator[j]:  # it means we're up or staying
                if switch_type[j] == 0:  # this means just go in the direction
                    running_down = False
                    bottom_length += switch_paths[j][int(running_down)]
                else:  # if switch_type is 1
                    # running_down stays the same
                    bottom_length += switch_paths[j][int(running_down)]
            elif running_state_bottom[0:num_values] == indicator_not[j]:  # if it doesn't match the indicator, as in it's down or switching
                if switch_type[j] == 0:  # this means it just goes down
                    running_down = True
                    bottom_length += switch_paths[j][int(running_down)]
                else:  # this means it's switching
                    running_down = not running_down
                    bottom_length += switch_paths[j][int(running_down)]
            else:
                print("ERROR: Check indicator array with state_table array; state mismatch.")
            running_state_bottom = running_state_bottom[num_values:]

        # Check to see if all values have been accounted for from state_table
        if len(running_state_bottom) != 0:
            # If this error is printed during runtime, something defined in this function is wrong
            print("ERROR: Not all states have been accounted for in bottom_path; state mismatch")

        # Calculate the path difference from the values accumulated by top_length & bottom_length
        path_diff = bottom_length - top_length
        # Can uncomment the line below to troubleshoot the bottom_ and top_lengths accumulated
        # print("path_diff: "+str(path_diff)+", top_length: "+str(top_length)+", bottom_length; "+str(bottom_length))
        # after all configs of state_table have been calculated, there should be a path difference for every config
        path_differences.append(path_diff)

    return np.array(path_differences)


# timing the full sweep for reference
big_start = time.time()

# create state_table and state_header for storage in pickle file
state_table, state_header = createStateTable(switch_table_file)

# calculate the path differences, one path difference for each configuration
path_lengths = calculatePathLengths(state_table, top_path, bottom_path, switch_paths)

# runs parameter sweep using switch info and returns results
all_results = runSweep(simulation_file, ONA_name, "dft_sweep", state_table, state_header, results_dict, top_path, bottom_path, switch_paths, min_wavelength, max_wavelength, num_segments, phase_shifts)

# store all of the results of the sweep in a pickle file for subsequent analysis
os.makedirs(os.path.dirname(pickle_file), exist_ok=True)  # make sure directory is created
with open(pickle_file, 'wb') as f:
    # order of variable storage and retrieval matters, so pickle filename can be very helpful for variable ordering
    pickle.dump([all_results, path_lengths, state_table, state_header], f)

# calculate total run time for reference
big_stop = time.time()
m, s = divmod(big_stop - big_start, 60)
h, m = divmod(m, 60)
print("Sweep and storage completed: "+str(h)+" hr, "+str(m)+" min.")
