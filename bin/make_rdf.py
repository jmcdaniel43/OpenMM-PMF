#!/usr/bin/env python2

from __future__ import print_function
import fileinput
import os
import sys
import re
from os import path
from subprocess import call
from sys import exc_info
import traceback
import argparse
from multiprocessing import Pool

from rdf import rdf
from rdf_smoother import smooth
from common import DiffusionSystem

startingDir = os.getcwd()

def make_rdf(system, rdfOutputFile, ion = False, bulk = False):
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

        if ion:
            solvent_resname = {
                "BF4": ("TME", "N"),
                "TMA": ("BF4", "B"),
                "TMEA": ("BF4", "B")
            }[diff_sys.diffusingIon]

            if diff_sys.ion_pair == "bmim":
                solvent_resname = ("BMI", "N1")
        else:
            solvent_resname = {
                "dce": ("dch", "CT"),
                "acn": ("acn", "CT")
            }[diff_sys.solvent]


        ion_atomselection = "(resname %s and name %s)" % (diff_sys.diffusingIon[:3], ion_atom)
        solvent_atomselection = "(resname %s and name %s)" % (solvent_resname[0][:3], solvent_resname[1])

        out = rdf(topology, trajectory, ion_atomselection, solvent_atomselection, bulk = bulk)

    except:
        print("error calculating rdf:", system)
        traceback.print_exc()
        return False

    with open(rdfOutputFile, "w") as f:
        for i in smooth(out):
            f.write(str(i) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ion", action="store_true", help="calculate ion coordination instead of solvation")
    parser.add_argument("--bulk", action="store_true", help="calculate bulk solvation coordination numbers")
    args = parser.parse_args()

    systems = []
    for line in sys.stdin.readlines():
        systems.append(line.strip())

    def make_ion_coordination(system):
        make_rdf(system, "ion_coordination.dat", ion = True)

    def make_bulk_solv(system):
        make_rdf(system, "bulk_solvation.dat", bulk = True)
    def make_bulk_ion(system):
        make_rdf(system, "bulk_ion.dat", ion = True, bulk = True)

    def __make_rdf(system):
        make_rdf(system, "rdf.dat")

    p = Pool(12)
    if args.ion:
        if args.bulk:
            p.map(make_bulk_ion, systems)
        else:
            p.map(make_ion_coordination, systems)
    else:
        if args.bulk:
            p.map(make_bulk_solv, systems)
        else:
            p.map(__make_rdf, systems)
