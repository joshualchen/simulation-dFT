import importlib.util

# This script will change all properties of the indicated component in all components of the simulation file,
# including within compound elements.
simulation_file = "simulation_dfts_basic.icp"
component_type = "WGD"
prop = "load from file"
setting = 0

# importing lumerical api
spec = importlib.util.spec_from_file_location("lumapi.py", "/Program Files/Lumerical/v211/api/python/lumapi.py")
lumapi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lumapi)
ic = lumapi.INTERCONNECT(filename=simulation_file, hide=True)
ic.switchtodesign()


def edit_node(base_element, component_type, prop, setting):
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
            edit_node(new_node, component_type, prop, setting)
            ic.groupscope(base_element)  # bring it back after completing inside compound element
            # compile the compound ones to check later in a string
        else:
            if ic.getnamed(element, "PREFIX") == component_type:
                ic.setnamed(element, prop, setting)


edit_node("::Root Element", component_type, prop, setting)
ic.save(simulation_file)
print("Completed")
