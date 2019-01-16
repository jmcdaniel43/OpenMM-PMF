#!/usr/bin/env python

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

def make_rdf(system, rdfOutputFile, ion = False, bulk = False, force = False):
    os.chdir(startingDir)

    if not path.exists(system):
        print("no system exists:", system)
        return []

    if path.exists(system + "/" + rdfOutputFile + ".dat"):
        print("coord exists:", system, rdfOutputFile)
        if force:
            print("overwriting")
        else:
            return []

    os.chdir(system)

    topology =  "md_nvt_prod_start_drudes.pdb"
    trajectory = "md_nvt_prod.dcd"

    diff_sys = DiffusionSystem(system)

    ion_atom = {
        "BF4": "B" ,
        "TMA": "N",
        "TMEA": "N"
    }[diff_sys.diffusingIon]

    if ion:
        solvent_resname = {
            "BF4": ("TME or resname TMA", "N"),
            "TMA": ("BF4", "B"),
            "TMEA": ("BF4", "B")
        }[diff_sys.diffusingIon]

        if diff_sys.ion_pair == "bmim":
            solvent_resname = ("BMI", "N1")
    else:
        solvent_resname = {
            "dce": ("dch", "CT CT1"),
            "acn": ("acn", "CT"),
            "h2o": ("HOH", "O"),
        }[diff_sys.solvent]


    ion_atomselection = "((resname %s) and name %s)" % (diff_sys.diffusingIon[:3], ion_atom)
    solvent_atomselection = "((resname %s) and name %s)" % (solvent_resname[0][:3], solvent_resname[1])

    coord = rdf(topology, trajectory, ion_atomselection, solvent_atomselection, bulk = bulk)

    smooth_coords = smooth(coord)

    with open(rdfOutputFile + ".dat", "w") as f:
        for i in smooth_coords:
            f.write(str(i) + "\n")

    return smooth_coords

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="force overwrite of exisiting files")
    args = parser.parse_args()
    print(args)

    systems = []
    for line in sys.stdin.readlines():
        systems.append(line.strip())

    def __rdf(system):
        make_rdf(system, "solv", force=args.force)
        make_rdf(system, "ion_coordination", force=args.force, ion=True)
        make_rdf(system, "bulk_solvation", force=args.force, bulk=True)
        make_rdf(system, "bulk_ion", force=args.force, bulk=True, ion=True)

    p = Pool(12)
    p.map(__rdf, systems)
