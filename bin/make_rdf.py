#!/usr/bin/env python

from __future__ import print_function
import fileinput
import os
import re
from os import path
from subprocess import call
from sys import exc_info
import traceback
from multiprocessing import Pool

from rdf import rdf
from rdf_smoother import smooth
from common import DiffusionSystem

startingDir = os.getcwd()

rdfOutputFile = "ion_coordination.dat"

def make_rdf(system):
    os.chdir(startingDir)

    if not path.exists(system):
        print("no system exists:", system)
        return False

    os.chdir(system)

    topology =  "md_nvt_prod_start_drudes.pdb"
    trajectory = "md_nvt_prod.dcd"

    diff_sys = DiffusionSystem(system)

    try:
        ion_atom = {
            "BF4": "B" ,
            "TMA": "N",
            "TMEA": "N"
        }[diff_sys.diffusingIon]

        solvent_resname = {
            "BF4": "TME",
            "TMA": "BF4",
            "TMEA": "BF4" 
        }[diff_sys.diffusingIon]

        ion_atomselection = "resname %s and name %s" % (diff_sys.diffusingIon[:3], ion_atom)
        solvent_atomselection = "(resname %s) and (name CT or name B or name N)" % (solvent_resname[:3])

        out = rdf(topology, trajectory, ion_atomselection, solvent_atomselection)

    except:
        print("error calculating rdf:", system)
        traceback.print_exc()
        return False

    with open(rdfOutputFile, "w") as f:
        for i in smooth(out):
            f.write(str(i) + "\n")

if __name__ == "__main__":
    p = Pool(12)
    p.map(make_rdf, map(lambda x: x.strip(), fileinput.input()))
