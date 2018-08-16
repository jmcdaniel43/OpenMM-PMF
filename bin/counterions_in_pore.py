#!/usr/bin/env python


from __future__ import division, print_function
from MDAnalysis import *
from MDAnalysis.lib.mdamath import triclinic_vectors
from MDAnalysis.analysis.distances import dist
from MDAnalysis.core.groups import AtomGroup

import numpy as np
from numpy import linalg

from multiprocessing import Pool
import argparse
import re

from minimum_image_displacement import minimum_image_disp

residRegex = re.compile("Ion index: \d+ Resid: (\d+)")

def count_counterions_in_pore(
        topology,
        trajectory,
        ion_atomselection,
        counterion_atomselection,
        pore_extent = 0,
        resid = "1",
        nframes = 6000,
        nwindows = 15):
    """
    Counts the number of counterions that are within the pore during specified windows.

    
    Parameters
    ----------

    pore_extent: the distance (in Angstroms) outside of the space between graphene sheets
        that is considered "inside" the pore
    """

    u = Universe(topology, trajectory)
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
    counterion_group = u.select_atoms(counterion_atomselection)
    pore_group = u.select_atoms("name C354 C274 C148 C68")

    # print(ion_group)
    # print(len(counterion_group))

    counts = []
    n_frames_per_window = nframes // nwindows
    for start in range(framestart, nframes, n_frames_per_window):
        c = count_counterions(u, ion_group, counterion_group, pore_group, start, n_frames_per_window)
        counts.append(c)

    return counts

def count_counterions(u, ion_group, counterion_group, pore_group, framestart, nframes_per_window):
    count = 0


    for i, frame in enumerate(u.trajectory[framestart : framestart + nframes_per_window]):
        box_vectors = np.mat(triclinic_vectors(u.trajectory.dimensions))
        box_inv = linalg.inv(box_vectors)

        reference_atom = pore_group[0].position
        z_start = reference_atom
        z_end = reference_atom
        z_start_atom = pore_group[0]
        z_end_atom = pore_group[0]

        for atom in pore_group[1:]:
            next_atom_minimum_disp = minimum_image_disp(box_vectors, reference_atom - atom.position, box_inv)
            next_atom_pos = reference_atom - next_atom_minimum_disp

            if next_atom_pos[2] < z_start[2]:
                z_start = next_atom_pos
                z_start_aotm = atom
            if next_atom_pos[2] > z_end[2]:
                z_end = next_atom_pos
                z_end_atom = atom

        z_center_of_pore = (z_start + z_end) / 2 # these two atoms' positions are in the same periodic image

        # print(z_start[2], z_end[2], z_center_of_pore[2])
        print_ion = False
        for cion in counterion_group:
            z_pos = reference_atom - minimum_image_disp(box_vectors, reference_atom - cion.position, box_inv)
            if z_pos[2] > z_start[2] and z_pos[2] < z_end[2]:
                # print(z_pos - z_start)
                count += 1
                print_ion = True

        if print_ion:
            ion_minimum_image_disp = minimum_image_disp(box_vectors, ion_group[0].position - z_center_of_pore, box_inv)
            print(ion_minimum_image_disp[2])

    return count

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("datadir", default="pmf_output", help="directory containing the output data from the simulation")
    parser.add_argument("ion_atomselection", default="resname BF4 and name B", help="atomselection string for the type of ion diffusing (resid 1 will be added to whatever is input)")
    parser.add_argument("counterion_atomselection", default="resname acn and name CT", help="atomselection string for the type of counterion")
    parser.add_argument("--ion_resid", default="1", help="ion residue id")
    args = parser.parse_args()

    # set the pdb topology and dcd trajectory
    topology = args.datadir + "/md_nvt_prod_start_drudes.pdb"
    trajectory = args.datadir + "/md_nvt_prod.dcd"

    count_counterions_in_pore(topology, trajectory, args.ion_atomselection, args.counterion_atomselection)
