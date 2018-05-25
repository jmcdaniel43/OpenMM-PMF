#!/usr/bin/env python

from __future__ import print_function
import fileinput
import os
import re
from os import path
from subprocess import call
from sys import exc_info
from traceback import print_exc

from rdf import rdf
from rdf_smoother import smooth
from common import DiffusionSystem

startingDir = os.getcwd()

rdfOutputFile = "rdf.log"

for sys in map(lambda x: x.strip(), fileinput.input()):
    os.chdir(startingDir)

    if not path.exists(sys):
        print("no system exists:", sys)
        continue

    os.chdir(sys)

    topology =  "md_nvt_prod_start_drudes.pdb"
    trajectory = "md_nvt_prod.dcd"

    diff_sys = DiffusionSystem(sys)

    try:
        ion_atom = {
            "BF4": "B" ,
            "TMA": "N",
            "TMEA": "N"
        }[diff_sys.diffusingIon]

        solvent_resname = {
            "acn": "acn",
            "dce": "dch" 
        }[diff_sys.solvent]

        ion_atomselection = "resname %s and name %s" % (diff_sys.diffusingIon[:3], ion_atom)
        solvent_atomselection = "resname %s and name CT" % (solvent_resname[:3])

        out = rdf(topology, trajectory, ion_atomselection, solvent_atomselection)

    except:
        print("error calculating rdf:", sys)
        print_exc()
        continue

    with open(rdfOutputFile, "w") as f:
        for i in smooth(out):
            f.write(str(i) + "\n")
