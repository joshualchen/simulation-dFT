import csv
import numpy as np
import pickle
import time
import os
import sys
sys.path.append("C:/Program Files/Lumerical/v212/api/python")  # may need to be changed for your personal directory
import lumapi

# Switch files below to save the output results from sweep
simulation_file = "simulation_dfts_basic_64.icp"  # icp file referenced by switch_table_file
switch_table_file = "switch_table_64.csv"  # csv file formatted according to guidelines
pickle_file = 'stored/803_64_2_[fsrmean_fsrmiddle_T]_pathlengths_statetable_stateheader.pickle'  # stored variables location
# Calculate path differences using top_path & bottom_path & switch_paths in function calculatePathLengths()
top_path = [0, 1, 2]  # these array elements refer to indices of switch_paths & other lists defined in calculatePathLengths()
bottom_path = [5, 6, 7]  # DEPENDING ON CIRCUIT DESIGN, PATH DIFFERENCE MAY BE CALCULATED DIFFERENTLY
switch_paths = np.array([[40, 0], [160, 0], [640, 0], [2560, 0], [10240, 0], [40, 120], [160, 480], [640, 1920],
                         [2560, 7680], [10240, 30720]])  # defines path lengths for each switch: [path up, path down]
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
        "Name": "T1",
        "Result": "::Root Element::ONA_3::input 1/mode 1/transmission",
        "Estimation": True,
        "Min": 1.0e11,
        "Max": 2.0e11
    },

    {
        "Name": "T2",
        "Result": "::Root Element::ONA_3::input 2/mode 1/transmission",
        "Estimation": True,
        "Min": 1.0e11,
        "Max": 2.0e11
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
            np_states = np.array([[float(s.strip()) for s in state.strip().strip("[]").split(",")] for state in states])
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
# CURRENTLY unused, setting numPoints manually in runSweep()
def numPoints(center_freq, freq_range, group_index, top_path, bottom_path, switch_paths):
    # center_freq/freq_range in Hz, switch_paths in um, theoretical FSR in m
    smallest_top_path = 0
    for switch in top_path:
        smallest_top_path += min(switch_paths[switch])
    largest_bottom_path = 0
    for switch in bottom_path:
        largest_bottom_path += max(switch_paths[switch])
    # calculating largest path difference; correlates to smallest FSR
    max_path_diff = largest_bottom_path - smallest_top_path  # paths in um
    # calculating the smallest wavelength given frequencies and speed of light
    min_wavelength = 2.998 * 10 ** 8 / (center_freq + freq_range / 2)  # wavelength in m
    # using the min wavelength because FSR to become a minimum
    theoretical_FSR = min_wavelength ** 2 / (group_index * max_path_diff) * 10 ** 6  # FSR in m
    # calculating the max wavelength to calculate the wavelength range
    max_wavelength = 2.998 * 10 ** 8 / (center_freq - freq_range / 2)
    # the necessary number of points to span the entire band, at least 10 points per FSR
    num_points = round((max_wavelength - min_wavelength) / theoretical_FSR * 10)
    return num_points


# opens Lumerical and executes the sweep, returns the results of the sweep
# inputs: .icp project file, the name of the sweep, the csv file with switch specifications, and the results to measure
def runSweep(project_file, sweep_name, switch_table_file, results_dict, top_path, bottom_path, switch_paths):
    # creating state_table and state_header using function createStateTable
    state_table, state_header = createStateTable(switch_table_file)

    ic = lumapi.INTERCONNECT(filename=project_file, hide=False)
    ic.switchtodesign()

    # adding FSR_mean & FSR_middle to the root element to sweep, must create new property
    ic.select("::Root Element")
    # string set to property "analysis script" is lumerical code to create the FSR_mean & FSR_middle properties
    ic.set("analysis script", 'FSR_data = getresult("::Root Element::ONA_3","input 1/mode 1/peak/free spectral range");'
                              + ' ?FSR_mean = mean(FSR_data.getattribute("mode 1 free spectral range (m)")); '
                              + 'setresult("FSR_mean", FSR_mean); '
                              + 'FSR_list = FSR_data.getattribute("mode 1 free spectral range (m)"); '
                              + '?FSR_middle = FSR_list(length(FSR_list)/2); setresult("FSR_middle", FSR_middle);')

    # change the number of points according to the simulation for efficiency, from the numPoints() function
    #center_freq = ic.getnamed("ONA_3", "center frequency")
    #freq_range = ic.getnamed("ONA_3", "frequency range")
    #group_index = 2  # May be subject to change, may need to move up for easy access as well
    #num_points = numPoints(center_freq, freq_range, group_index, top_path, bottom_path, switch_paths)
    #ic.setnamed("ONA_3", "number of points", num_points)

    # measure the amount of time for one run to predict total run time
    start = time.time()
    ic.run()
    end = time.time()
    run_time_one = end - start

    # setup sweep
    # delete previous instances of sweep
    ic.deletesweep(sweep_name)
    ic.addsweep(0)
    ic.setsweep("sweep", "name", sweep_name)
    ic.setsweep(sweep_name, "type", "Values")
    ic.setsweep(sweep_name, "number of points", state_table.shape[0])
    print("Sweep created: " + ic.getsweep(sweep_name, "name"))

    # setup sweep parameters according to state_table and state_header
    # add all sweep parameters as dictionaries
    for i in range(len(state_header)):
        object_name = state_header[i][0]
        property_name = state_header[i][1]
        para = {
            "Name": object_name + ", " + property_name,
            "Parameter": "::Root Element::" + object_name + "::" + property_name,
            "Type": "Number"
        }
        # fill the parameter settings for each object/property combo using data from state_table
        for j in range(state_table.shape[0]):
            para["Value_"+str(j+1)] = np.transpose(state_table)[i][j]
        # add the created dictionary as another sweep parameter
        ic.addsweepparameter(sweep_name, para)

    # set up results parameters
    for result in results_dict:
        ic.addsweepresult(sweep_name, result)

    run_time_total = state_table.shape[0] * run_time_one
    m, s = divmod(run_time_total, 60)
    h, m = divmod(m, 60)
    print("Running sweep. Crude time estimate: "+str(h)+" hr, "+str(m)+" min")
    ic.runsweep(sweep_name)
    print("Sweep completed.")

    # get sweep results and return in a list
    results = []
    for result in results_dict:
        results.append(ic.getsweepresult(sweep_name, result["Name"]))

    return results


# DEFINING A FUNCTION TO RUN ONE SWEEP based on the voltage setting
# This function is helpful to troubleshoot a specific state (row of state_table) if sweep results are not as expected
def runOne(project_file, switch_state, state_header):
    # opening lumerical
    spec = importlib.util.spec_from_file_location("lumapi.py", "/Program Files/Lumerical/v211/api/python/lumapi.py")
    lumapi = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(lumapi)
    ic = lumapi.INTERCONNECT(filename=project_file, hide=False)
    ic.switchtodesign()
    # setting the voltage states in the simulation according to switch_state
    for i in range(len(state_header)):
        object_name = state_header[i][0]
        property_name = state_header[i][1]
        ic.setnamed(object_name, property_name, switch_state[i])
    # running the simulation
    ic.run()


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
    switch_type = [0, 1, 1, 1, 1, 0, 1, 1, 1, 1]
    # indicator also mirrors indices of switch_paths & switch_type
    # indicator lists for each switch the setting for "go up" for switch_type of 0, or "stay" for switch_type of 1; i.e. "State 1"
    indicator = [[-2.35,-0.5],[-0.5,-3.593],[-0.5,-3.593],[-0.5,-3.593],[-0.5,-3.593],[-2.35,-0.5],[-0.5,-3.593],[-0.5,-3.593],[-0.5,-3.593],[-0.5,-3.593]] # up or stay indicator
    # indicator_not lists the setting for "go down" for switch_type of 0, or "switch" for switch_type of 1; i.e. "State 2"
    indicator_not = [[-0.5,-2.35],[-0.5,-0.5],[-0.5,-0.5],[-0.5,-0.5],[-0.5,-0.5],[-0.5,-2.35],[-0.5,-0.5],[-0.5,-0.5],[-0.5,-0.5],[-0.5,-0.5]] # down or switch indicator

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
            if running_state_top[0:num_values] == indicator[j]: # it means we're up or staying
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


if __name__ == "__main__":
    # create state_table and state_header for storage in pickle file
    state_table, state_header = createStateTable(switch_table_file)

    # calculate the path differences, one path difference for each configuration
    path_lengths = calculatePathLengths(state_table, top_path, bottom_path, switch_paths)

    # runs parameter sweep using switch info and returns results
    results = runSweep(simulation_file, "run_dft_sweep", switch_table_file, results_dict, top_path, bottom_path, switch_paths)

    # store all of the results of the sweep in a pickle file for subsequent analysis
    os.makedirs(os.path.dirname(pickle_file), exist_ok=True)  # make sure directory is created
    with open(pickle_file, 'wb') as f:
        # order of variable storage and retrieval matters, so pickle filename can be very helpful for variable ordering
        pickle.dump([results, path_lengths, state_table, state_header], f)
