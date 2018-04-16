#!/usr/bin/env python

# this script uses the MD Analysis tools to compute RDFs.  MD Analysis should be loaded with following lines
#    module load anaconda3/latest
#    source activate p4env
#

from __future__ import division, print_function # python2 defaults to interger division. py3 defaults to floating division
from MDAnalysis import *
import MDAnalysis.analysis.rdf as rdf
from MDAnalysis.lib.mdamath import triclinic_vectors, box_volume
import numpy as np
import scipy as sp
import argparse
from math import cos

parser = argparse.ArgumentParser()
#parser.add_argument("datadir", help="path where the {1sshet,2sheet} directories are")
parser.add_argument("datadir", default="pmf_output", help="directory containing the output data from the simulation")
parser.add_argument("--rdf_shell", default=20 ,type=int, help="size (in Angstroms) of solvation shell over which to calculate RDF")
parser.add_argument("--is", "--integration_shell", default=6, type=int, help="size (in Angstroms) of solvation shell over which to integrate", dest="integration_shell")
parser.add_argument("--bulk", default="false", help="whether to calculate the bulk rdf values")
args = parser.parse_args()

# set the atoms for computing RDF
atoms1=('B',)

# set the pdb topology and dcd trajectory
topology = args.datadir + "/md_nvt_prod_start_drudes.pdb"
trajectory = args.datadir + "/md_nvt_prod.dcd"

u=Universe(topology, trajectory)

# frame to start
framestart=u.trajectory.n_frames - 1199

group1 = u.select_atoms("resname BF4 and resid 1 and name B")
group2 = u.select_atoms("name CT and resname acn")

if args.bulk.lower() in ['t', 'tr', 'tru', 'true']:
    for ele in atoms1:
        group1 = group1 + u.select_atoms("name %s" % ele)

def makeRDF(startFrame):
    frame = startFrame + framestart
    rdff = rdf.InterRDF(group1 , group2 , nbins=100, range=(0.0, args.rdf_shell), start=frame, stop=frame+19)
    rdff.run()
    return rdff

rdfs = list(map(makeRDF, range(0,1200,20)))

end = int(len(rdfs[0].bins) * (args.integration_shell / args.rdf_shell))

dims = u.trajectory.dimensions
#boxVecs = triclinic_vectors(dims)
#boxVolume = np.dot(np.cross(boxVecs[0], boxVecs[1]), boxVecs[2])
boxVolume = box_volume(dims)
density = group2.n_atoms / boxVolume

for rdf in rdfs:
    data = [g * rdf.bins[i]**2 for i, g in enumerate(rdf.rdf[:end])]
    coordinationNum = (4 * 3.14159 * density) * sp.integrate.trapz(data, rdf.bins[:len(data)])
    print(coordinationNum)
    # for i in range(len(rdf.bins)):
    #     print("%f %f" % (rdf.bins[i], rdf.rdf[i]), end="")
    # print()

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
