#!/usr/bin/env python

from __future__ import division, print_function # python2 defaults to interger division. py3 defaults to floating division
from MDAnalysis import *
from MDAnalysis.analysis import rdf as mdaRDF
from MDAnalysis.lib.mdamath import triclinic_vectors, box_volume
import numpy as np
import scipy as sp
import argparse
import re

from rdf_smoother import smooth

residRegex = re.compile("Ion index: \d+ Resid: (\d+)")

def rdf(topology,
        trajectory,
        ion_atomselection,
        solvent_atomselection,
        rdf_shell = 20,
        integration_shell = 6,
        bulk = False,
        print_rdf = None,
        no_smooth = False,
        resid = "1",
        nframes = 6000):

    u=Universe(topology, trajectory)
    framestart = u.trajectory.n_frames - nframes

    # for the newer simulations, the resid is no longer always 1
    # so we can check for the resid in those new sims' logs
    # if the regex matches something, it's a new sim: record the resid
    # otherwise, keep using 1 as the resid
    with open("output.log") as log:
        i = 0
        for line in log.readlines():
            if i > 200:
                break # there probably is no match for the regex

            m = residRegex.match(line)
            if m is not None:
                resid = m.group(1)
                break

            i += 1

    ion_group = u.select_atoms(ion_atomselection + " and resid " + resid)
    solvent_group = u.select_atoms(solvent_atomselection)

    # print(ion_group)
    # print(solvent_group)

    if bulk:
        ion_group += u.select_atoms(ion_atomselection)

    rdfs = []
    window_frames = nframes // 60
    for i in range(0, nframes, window_frames):
        frame = i + framestart
        rdf = mdaRDF.InterRDF(ion_group, solvent_group, nbins=100, range=(0.0, rdf_shell), start=frame, stop=frame + window_frames)
        rdf.run()
        rdfs.append(rdf)

    end = int(len(rdfs[0].bins) * (integration_shell / rdf_shell))

    dims = u.trajectory.dimensions
    boxVolume = box_volume(dims)
    density = solvent_group.n_atoms / boxVolume

    # print(boxVolume, density)

    return_strings = []

    if print_rdf is not None:
        for i in range(len(rdfs[print_rdf].bins)):
            return_strings.append("%f %f" % (rdfs[print_rdf].bins[i], rdfs[print_rdf].rdf[i]))

    else:
        for rdf in rdfs:
            data = [g * rdf.bins[i]**2 for i, g in enumerate(rdf.rdf[:end])]
            coordinationNum = (4 * 3.14159 * density) * sp.integrate.trapz(data, rdf.bins[:len(data)])
            return_strings.append(coordinationNum)
    
    if no_smooth or (print_rdf is not None):
        return return_strings

    return smooth(map(float, return_strings))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("datadir", default="pmf_output", help="directory containing the output data from the simulation")
    parser.add_argument("ion_atomselection", default="resname BF4 and name B", help="atomselection string for the type of ion diffusing (resid 1 will be added to whatever is input)")
    parser.add_argument("solvent_atomselection", default="resname acn and name CT", help="atomselection string for the type of solvent")
    parser.add_argument("--rdf_shell", default=20 ,type=int, help="size (in Angstroms) of solvation shell over which to calculate RDF")
    parser.add_argument("--is", "--integration_shell", default=6, type=int, help="size (in Angstroms) of solvation shell over which to integrate", dest="integration_shell")
    parser.add_argument("--bulk", action="store_true", help="whether to calculate the bulk rdf values")
    parser.add_argument("--print_rdf", default=None, type=int, help="print rdf values for frame <x>; default print coordination number")
    parser.add_argument("--no_smooth", action="store_true", help="do not smooth coordination number output")
    args = parser.parse_args()

    # set the pdb topology and dcd trajectory
    topology = args.datadir + "/md_nvt_prod_start_drudes.pdb"
    trajectory = args.datadir + "/md_nvt_prod.dcd"

    out = rdf(topology, trajectory, args.ion_atomselection, args.solvent_atomselection,
            bulk=args.bulk,
            integration_shell=args.integration_shell,
            rdf_shell=args.rdf_shell,
            print_rdf=args.print_rdf,
            no_smooth=args.no_smooth)

    for i in out:
        print(i)

# MANUAL CALCULATION
# coordNums = [0] * 60
# for i, dz in enumerate(range(0,1200,20)):
#     for j in range(20):
#         u.trajectory[20 * i + j]
#         coordNum = 0
#         group2 = u.select_atoms("(around 6 (name B and resid 1)) and (name CT)")
#         for atom in group2:
#             coordNum += 1 / group1.n_atoms
#         coordNums[i] += coordNum / 20
# 
# for i in coordNums:
#     print(i)
