#!/usr/bin/env python

import re
import sys
from pathlib import Path
from math import sqrt, floor, ceil
import multiprocessing as mp
from multiprocessing import Pool
import argparse
import traceback

import pandas as pd
import numpy as np
from MDAnalysis import *

from common import DiffusionSystem
from make_pmf import make_pmf
from rdf import make_rdf
from calc_dipole import calc_dipole
from convert_umbrella_output import *

"""
Retrieves data from PMF and RDF files
"""


graph_atoms = {
    10: ['C38', 'C48', 'C58', 'C139', 'C149', 'C159', 'C169', 'C230', 'C231', 'C239', 'C242', 'C253', 'C263', 'C269', 'C270', 'C273', 'C281', 'C282', 'C283', 'C330', 'C331', 'C342', 'C353', 'C364', 'C370', 'C374', 'C381', 'C382', 'C383', 'C384'],
    14: ['C27', 'C28', 'C37', 'C47', 'C57', 'C68', 'C128', 'C129', 'C138', 'C148', 'C158', 'C168', 'C179', 'C220', 'C221', 'C229', 'C232', 'C243', 'C254', 'C264', 'C274', 'C279', 'C280', 'C284', 'C320', 'C321', 'C332', 'C343', 'C354', 'C365', 'C375', 'C380', 'C385', 'C691', 'C692', 'C693', 'C694', 'C791', 'C792', 'C793', 'C794', 'C795']
}


pmf_regex = re.compile("(-?\d+\.\d+)\s+(\d+\.\d+)\s+.*")

class PMF_Data:
    def __init__(self, data):
        self.start = list(data.values())[0]
        self.max = max(data.values())
        self.end = list(data.values())[-1]
        self.data = data

    def __len__(self):
        return len(self.data)

class RDF_Data:
    def __init__(self, data):
        self.start = data[0]
        self.min = min(data)
        self.end = data[-1]
        self.data = data
        self.average = sum(data) / len(data)

    def __len__(self):
        return len(self.data)

class DiffusionData:
    def __init__(self, pandas_frame, force = False):
        self.pmf_files = {}
        self.solv_files = {}
        self.ion_files = {}
        self.bulk_solv_files = {}
        self.bulk_ion_files = {}
        self.graph = {}

        try:
            if force:
                raise FileNotFoundError("force recreation")

            self.df = pd.read_hdf(pandas_frame, 'porous')
        except FileNotFoundError:
            row_index = pd.MultiIndex.from_product([[1,2], [7,10,14], ['acn','dce'], ['bmim','tma_tmea'], ['BF4', 'TMA', 'TMEA']])
            self.df = pd.DataFrame(index=row_index, columns=['graph', 'z', 'pmf', 'solv', 'ion'])
            
            idx = pd.IndexSlice
            # drop rows that are invalid
            for sheet in [1,2]:
                for solv in ["acn", "dce"]:
                    for ion in ["TMA", "TMEA"]:
                        self.df = self.df.drop(idx[(sheet, 7, solv, "tma_tmea", ion)])
                        for pore in [7,10,14]:
                            self.df = self.df.drop(idx[(sheet, pore, solv, "bmim", ion)])

    def check_present(key):
        column = self.df.loc(key)
        for i in column:
            pass

        if key in self.df.loc(key):
            print("data present", system)
            return True

        return False

    def apply_file(self, path, force = False):
        idx = pd.IndexSlice
        system = DiffusionSystem(path)
        key = idx[system.pandas_tuple()]

        try:
            if np.isnan(self.df.loc[key].graph).any():
                u = Universe(path+"/md_nvt_prod_start_drudes.pdb", path+"/md_nvt_prod.dcd")
                with open(path+"/output.log") as f:
                    pos = electrode_positions(f.readlines(), u)
                    print(system)
                    print(pos)
                    self.df.loc[key].graph = pos

            pmf = get_pmf(path)
            self.df.loc[key].z = list(pmf.data.keys())
            self.df.loc[key].pmf = list(pmf.data.values())

            # get all the rdf values
            rdfs = get_rdfs(system)
            self.df.loc[key].solv = rdfs["solv"]
            self.df.loc[key].ion = rdfs["ion"]
        except Exception as e:
            print("failed", system)
            traceback.print_exc()
            if type(e) == KeyboardInterrupt:
                return

def read_pmf_file(path):
    pmf_values = {}
    with open(path) as pmf_file:
        for line in pmf_file:
            pmf_value = pmf_regex.match(line)
            if pmf_value is not None:
                pmf_values[float(pmf_value.group(1))] = float(pmf_value.group(2))

    if len(pmf_values) < 1:
        return

    return PMF_Data(pmf_values)

def get_pmf(datapath):
    pmf_outputfile_path = Path(datapath) / "pmf" / "pmf"
    if not pmf_outputfile_path.exists():
        print("no such pmf file", pmf_path, file=sys.stderr)
        print("attempting to calculate", file=sys.stderr)
        if make_pmf(datapath) is not None:
            raise ValueError("pmf calculation failed")

    return read_pmf_file(pmf_outputfile_path)

def get_rdf(system, extension, ion = False, bulk = False):
    rdf_path = Path(system.sim_path() + "/" + extension + ".dat")
    rdf_values = []

    if rdf_path.exists():
        with open(rdf_path) as rdf_file:
            for line in rdf_file:
                rdf_values.append(float(line.strip()))

    if len(rdf_values) < 1:
        print("no such rdf file", rdf_path, file=sys.stderr)
        print("attempting to calculate", file=sys.stderr)
        rdf_values = make_rdf(system.sim_path(), extension, ion = ion, bulk = bulk)
        if len(rdf_values) < 1:
            raise ValueError("rdf calculation failed")

    return RDF_Data(rdf_values)

def get_rdfs(system):
    solv      = get_rdf(system, "solv", False, True)
    ion       = get_rdf(system, "ion", True, True)
    bulk_solv = get_rdf(system, "bulk_solv", False, False)
    bulk_ion  = get_rdf(system, "bulk_ion", True, False)

    return {"solv": (bulk_solv.average, solv.data), "ion": (bulk_ion.average, ion.data)}

def electrode_positions(output_lines, u):
    lines = output_lines
    for line in range(len(lines)): # find start and maximum z dimension from output file
        startingDistance_match = startingDistanceRegex.match(lines[line])
        startZ_match1 = startingPotentialDimensionRegex1.match(lines[line])
        startZ_match2 = startingPotentialDimensionRegex2.match(lines[line])
        startZ_match3 = startingPotentialDimensionRegex3.match(lines[line])

        if startingDistance_match is not None:
            startingDistance = float(startingDistance_match.group(1)) * 10

        elif startZ_match1 is not None:
            start = float(startZ_match1.group(5)) * 10 # convert nm to angstrom
        elif startZ_match2 is not None:
            start = float(startZ_match2.group(5)) * 10 # convert nm to angstrom
        elif startZ_match3 is not None:
            start = float(startZ_match3.group(5)) * 10 # convert nm to angstrom

        elif startRegex.match(lines[line]) is not None:
            line += 1
            break

    try:
        center_of_pore = start + startingDistance
        print("electrode_calc: startingDistance", startingDistance, "Å")
    except NameError:
        center_of_pore = start + 10 # usually simulations start 1 nm from the center of the pore
        print("electrode_calc: startingDistance", 10, "Å")

    g = []
    g_pos = [[], []]
    for i in range(2):
        # get the atoms in each residue of the electrode
        g.append(u.select_atoms("resid " + str(i+1) + " and resname grp"))

    for frame in u.trajectory:
        graph_positions = np.array(list(map(lambda x: x.positions,g)))
        mean_pos = np.mean(graph_positions, axis=1)[:,2]
        g_pos[0].append(mean_pos[0] - center_of_pore)
        g_pos[1].append(mean_pos[1] - center_of_pore)

    return np.mean(g_pos, axis=1).tolist()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", type=int, default=1, help="number of processes to launch to calculate things concurrently")
    parser.add_argument("--label", default='porous')
    parser.add_argument("--force_df", action='store_true')
    parser.add_argument("--force_files", action='store_true')
    args = parser.parse_args()

    filename = "data.hd5"
    dd = DiffusionData(filename, args.force_df)

    procs = []
    with Pool(processes=args.j) as pool:
        for i in sys.stdin.readlines():
            procs.append(pool.apply_async(dd.apply_file, [i.strip(), args.force_files]))
        for p in procs:
            p.get()

    dd.df.to_hdf(filename, "porous")
    print(dd.df)
