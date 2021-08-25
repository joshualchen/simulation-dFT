import importlib.util

# This script will change all "ldf filename" and "s parameters filename" in all components of the simulation file,
# including within compound elements. Also, running this script assumes all "load from file" options to become true,
# use change_property_all to change all back
simulation_file = "simulation_dfts_basic_256.icp"
old_dir = "C:/Users/linda/Downloads/For Dylan/components/"
new_dir = "C:/Users/linda/Desktop/For Josh/components/"

# importing lumerical api
spec = importlib.util.spec_from_file_location("lumapi.py", "/Program Files/Lumerical/v211/api/python/lumapi.py")
lumapi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lumapi)
ic = lumapi.INTERCONNECT(filename=simulation_file, hide=True)
ic.switchtodesign()


def edit_node(base_element, old_dir, new_dir):
    # setting to the input base_element
    ic.groupscope(base_element)
    components = ic.getdata()
    a = components.split("\n")
    # for each element in base_element
    for i in range(0, len(a)):
        # what is the name of the element? below
        element = a[i]
        if ic.getnamed(element, "PREFIX") == "COMPOUND":
            # go into the new compound element node
            new_node = base_element + "::" + element
            edit_node(new_node, old_dir, new_dir)
            ic.groupscope(base_element)  # bring it back after completing inside compound element
            # compile the compound ones to check later in a string
        else:
            # first, try to see if they have ldf filename
            try:
                ic.setnamed(element, "load from file", 1)
                old_ref = ic.getnamed(element, "ldf filename")  # renames the correct folder for ldf
                new_ref = old_ref.replace(old_dir, new_dir)
                ic.setnamed(element, "ldf filename", new_ref)
            except:
                continue
            # next, try to see if they have s parameters filename
            try:
                ic.setnamed(element, "load from file", 1)
                old_ref = ic.getnamed(element, "s parameters filename")
                new_ref = old_ref.replace(old_dir, new_dir)
                ic.setnamed(element, "s parameters filename", new_ref)
            except:
                continue


edit_node("::Root Element", old_dir, new_dir)
ic.save(simulation_file)
print("Completed")
