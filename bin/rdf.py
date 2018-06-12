#!/usr/bin/env python

from __future__ import division, print_function # python2 defaults to interger division. py3 defaults to floating division
from MDAnalysis import *
from MDAnalysis.analysis import rdf as mdaRDF
from MDAnalysis.lib.mdamath import triclinic_vectors, box_volume
import numpy as np
import scipy as sp
import argparse
from rdf_smoother import smooth

def rdf(topology,
        trajectory,
        ion_atomselection,
        solvent_atomselection,
        rdf_shell = 20,
        integration_shell = 6,
        bulk = False,
        print_rdf = False,
        no_smooth = False):

    u=Universe(topology, trajectory)
    framestart = u.trajectory.n_frames - 1199

    ion_group = u.select_atoms(ion_atomselection + " and resid 1")
    solvent_group = u.select_atoms(solvent_atomselection)

    # print(ion_atomselection)
    # print(solvent_atomselection)

    if bulk:
        ion_group += u.select_atoms(ion_atomselection)

    rdfs = []
    for i in range(0,1200,20):
        frame = i + framestart
        rdf = mdaRDF.InterRDF(ion_group , solvent_group , nbins=100, range=(0.0, rdf_shell), start=frame, stop=frame+19)
        rdf.run()
        rdfs.append(rdf)

    end = int(len(rdfs[0].bins) * (integration_shell / rdf_shell))

    dims = u.trajectory.dimensions
    boxVolume = box_volume(dims)
    density = solvent_group.n_atoms / boxVolume

    # print(boxVolume, density)

    return_strings = []

    for rdf in rdfs:
        if print_rdf:
            for i in range(len(rdf.bins)):
                return_strings.append("%f %f" % (rdf.bins[i], rdf.rdf[i]), end="")
        else:
            data = [g * rdf.bins[i]**2 for i, g in enumerate(rdf.rdf[:end])]
            coordinationNum = (4 * 3.14159 * density) * sp.integrate.trapz(data, rdf.bins[:len(data)])
            return_strings.append(coordinationNum)
    
    if no_smooth:
        return return_strings
    else:
        return smooth(map(float, return_strings))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("datadir", default="pmf_output", help="directory containing the output data from the simulation")
    parser.add_argument("ion_atomselection", default="resname BF4 and name B", help="atomselection string for the type of ion diffusing (resid 1 will be added to whatever is input)")
    parser.add_argument("solvent_atomselection", default="resname acn and name CT", help="atomselection string for the type of solvent")
    parser.add_argument("--rdf_shell", default=20 ,type=int, help="size (in Angstroms) of solvation shell over which to calculate RDF")
    parser.add_argument("--is", "--integration_shell", default=6, type=int, help="size (in Angstroms) of solvation shell over which to integrate", dest="integration_shell")
    parser.add_argument("--bulk", action="store_true", help="whether to calculate the bulk rdf values")
    parser.add_argument("--print_rdf", action="store_true", help="print rdf values; default print coordination number")
    parser.add_argument("--no_smooth", action="store_true", help="do not smooth coordination number output")
    args = parser.parse_args()

    # set the pdb topology and dcd trajectory
    topology = args.datadir + "/md_nvt_prod_start_drudes.pdb"
    trajectory = args.datadir + "/md_nvt_prod.dcd"

    out = rdf(topology, trajectory, **args)

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
