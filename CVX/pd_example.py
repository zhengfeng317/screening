#!/usr/bin/env python3

"""
This is a basic example of how to create, plot, and analyze Phase Diagrams using the pymatgen
codebase and Materials Project database. To run this example, you should:
* have pymatgen (www.pymatgen.org) installed along with matplotlib
* obtain a Materials Project API key (https://www.materialsproject.org/open)
* paste that API key in the MAPI_KEY variable below, e.g. MAPI_KEY = "foobar1234"
For citation, see https://www.materialsproject.org/citing
For the accompanying comic book, see http://www.hackingmaterials.com/pdcomic
"""

from pymatgen.ext.matproj import MPRester
from pymatgen.entries.compatibility import MaterialsProjectCompatibility
#from pymatgen.phasediagram.pdanalyzer import PDAnalyzer
#from pymatgen.phasediagram.pdmaker import PhaseDiagram
#from pymatgen.phasediagram.plotter import PDPlotter
#from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from pymatgen.analysis.phase_diagram import *

import sys

if __name__ == "__main__":
    MAPI_KEY = 'oW693Gka54dOv8rc'  # You must change this to your Materials API key! (or set MAPI_KEY env variable)
    #system = ["Ce", "B", "C"]  # system we want to get PD for
    system = [sys.argv[1], sys.argv[2], sys.argv[3]]

    mpr = MPRester(MAPI_KEY)  # object for connecting to MP Rest interface
    compat = MaterialsProjectCompatibility()  # sets energy corrections and +U/pseudopotential choice

    # Create phase diagram!
    unprocessed_entries = mpr.get_entries_in_chemsys(system)
    processed_entries = compat.process_entries(unprocessed_entries)  # filter and add energy corrections
    pd = PhaseDiagram(processed_entries)

    # Plot!
    #plotter = PDPlotter(pd, show_unstable=True)  # you can also try show_unstable=True
    #plotter.show()
    #plotter.write_image("{}.png".format('-'.join(system)), "png")  # save figure

    # Analyze phase diagram!
    #pda = PDAnalyzer(pd)

    print('Stable Entries (formula, materials_id)\n--------')
    #print(pd.stable_entries)
    with open("mp_stable_ene.dat", "w+") as fout:
        for e in pd.stable_entries:
            #print(e)
            print(e.composition.reduced_formula, e.entry_id, e.energy_per_atom)
            fout.write("%10s%15s%15.9f" %(e.composition.reduced_formula, e.entry_id, e.energy_per_atom))
            fout.write("\n")

    #print('\nUnstable Entries (formula, materials_id, e_above_hull (eV/atom), decomposes_to)\n--------')
    with open("mp_unstable_ene.dat", "w+") as fout2:
        for e in pd.unstable_entries:
            decomp, e_above_hull = pd.get_decomp_and_e_above_hull(e)
            pretty_decomp = [("{}:{}".format(k.composition.reduced_formula, k.entry_id), round(v, 2)) for k, v in decomp.items()]
            #print(e.composition.reduced_formula, e.entry_id, "%.3f" % e_above_hull, pretty_decomp)
            fout2.write("%10s%15s%15.9f" %(e.composition.reduced_formula, e.entry_id, e.energy_per_atom))
            fout2.write("\n")
