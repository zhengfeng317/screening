#!/usr/bin/env python3

import os
import shutil
import numpy as np
from pymatgen.ext.matproj import MPRester
#from matminer.featurizers.site import CrystalNNFingerprint
#from matminer.featurizers.structure import SiteStatsFingerprint
from pymatgen.core import Lattice, Structure, Molecule

MAPI_KEY = 'oW693Gka54dOv8rc'
mpr = MPRester(MAPI_KEY)

with open("./mp_stable_ene.dat", "r") as fin:
    lines = fin.readlines()
    for i in range(len(lines)):
        ll = lines[i].split()
        struc_com = ll[0]
        struc_id = ll[1]
        struc_mp = mpr.get_structure_by_material_id(struc_id)
        struc_mp.to(filename="POSCAR")
        pos_mp = "POSCAR." + struc_id + "." + struc_com
        shutil.copy("POSCAR", pos_mp)
        #os.system("mv POSCAR " + pos_mp)
