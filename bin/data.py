#!/usr/bin/env python

import re
import sys
from pathlib import Path
from math import sqrt, floor, ceil
from multiprocessing import Pool
import argparse

import pandas as pd
import numpy as np

from common import DiffusionSystem
from make_pmf import make_pmf
from make_rdf import make_rdf
from make_grace import graceScript
from calc_dipole import calc_dipole

"""
Retrieves data from PMF and RDF files
"""

pmf_regex = re.compile("(-?\d+\.\d+)\s+(\d+\.\d+)\s+.*")

def get_pmf(datapath):
    pmf_path = Path(datapath) / "pmf" / "pmf"
    if not pmf_path.exists():
        print("no such pmf file", datapath, file=sys.stderr)
        return
        # print("attempting to calculate", file=sys.stderr)
        # if make_pmf(datapath) is not None:
        #     print("pmf calculation failed", file=sys.stderr)
        #     return

    pmf_values = {}

    with open(pmf_path) as pmf_file:
        for line in pmf_file:
            pmf_value = pmf_regex.match(line)
            if pmf_value is not None:
                pmf_values[float(pmf_value.group(1))] = float(pmf_value.group(2))

    if len(pmf_values) < 1:
        return

    return PMF_Data(pmf_values)

def get_rdf(datapath):
    rdf_path = Path(datapath)
    if not rdf_path.exists():
        print("no such rdf file", rdf_path, file=sys.stderr)
        return
        # print("attempting to calculate", file=sys.stderr)
        # if make_rdf(datapath) is not None:
        #     print("pmf calculation failed", file=sys.stderr)
        #     return

    rdf_values = []

    with open(rdf_path) as rdf_file:
        for line in rdf_file:
            rdf_values.append(float(line.strip()))

    if len(rdf_values) < 1:
        return

    return RDF_Data(rdf_values)

def find_sys(inp, systems):
    try:
        system_filter = re.compile(inp)
    except:
        # probably an invalid regex
        print("regex excceptiong", file=sys.stderr)

    return [x for x in systems if system_filter.match(x.output_name()) is not None]

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
    def __init__(self, datapaths):
        self.pmf_files = {}
        self.solv_files = {}
        self.ion_files = {}
        self.bulk_solv_files = {}
        self.bulk_ion_files = {}

        for line in datapaths:
            path = line.strip()
            system = DiffusionSystem(path)

            pmf = get_pmf(path)
            if pmf is not None:
                for k, v in self.pmf_files.items():
                    if len(v) != len(pmf):
                        print("this pmf file is not like the others", path, "(compared to", k, ")", file=sys.stderr)
                    break

                self.pmf_files[system] = pmf

            solv = get_rdf(path + "/solv.dat")
            if solv is not None:
                for k, v in self.solv_files.items():
                    if len(v) != len(solv):
                        print("this rdf file is not like the others", path, "(compared to", k, ")", file=sys.stderr)
                    break

                self.solv_files[system] = solv

            ion = get_rdf(path + "/ion_coordination.dat")
            if ion is not None:
                for k, v in self.ion_files.items():
                    if len(v) != len(ion):
                        print("this ion file is not like the others", path, "(compared to", k, ")", file=sys.stderr)
                    break

                self.ion_files[system] = ion

            bulk_solv = get_rdf(path + "/bulk_solvation.dat")
            if bulk_solv is not None:
                for k, v in self.bulk_solv_files.items():
                    if len(v) != len(bulk_solv):
                        print("this bulk file is not like the others", path, "(compared to", k, ")", file=sys.stderr)
                    break

                self.bulk_solv_files[system] = bulk_solv

            bulk_ion = get_rdf(path + "/bulk_ion.dat")
            if bulk_ion is not None:
                for k, v in self.bulk_ion_files.items():
                    if len(v) != len(bulk_ion):
                        print("this bulk file is not like the others", path, "(compared to", k, ")", file=sys.stderr)
                    break

                self.bulk_ion_files[system] = bulk_ion

    def make_grace(self, systems):
        for system in systems:
            try:
                pmf = self.pmf_files[system].data
                coordinates = list(pmf.keys())
                energies = list(pmf.values())

                energiesStr = ""
                for z, e in pmf.items():
                    energiesStr += "{:f} {:f}\n".format(z, e)

                start_x = floor(coordinates[0])
                end_x = ceil(coordinates[-1])
                end_y = 20
                grace = graceScript.format(start_x, end_x, end_y, "pmf", "{:s}") + energiesStr + "&"

                with open("grace/pmf_" + system.output_name() + ".xgr", "w") as f:
                    f.write(grace.format(system.output_name()))
            except:
                print("pmf grace data failed")
                print(sys.exc_info())
                continue

            try:
                solv = self.solv_files[system].data
                energiesStr = ""
                for z, e in zip(coordinates, solv):
                    energiesStr += "{:f} {:f}\n".format(z, e)

                grace = graceScript.format(start_x, end_x, end_y, "solv", "{:s}") + energiesStr + "&"

                with open("grace/solv_" + system.output_name() + ".xgr", "w") as f:
                    f.write(grace.format(system.output_name()))
            except:
                print("solv grace data failed")
                print(sys.exc_info())

            try:
                ion = self.ion_files[system].data
                energiesStr = ""
                for z, e in zip(coordinates, ion):
                    energiesStr += "{:f} {:f}\n".format(z, e)

                grace = graceScript.format(start_x, end_x, end_y, "ion", "{:s}") + energiesStr + "&"

                with open("grace/ion_" + system.output_name() + ".xgr", "w") as f:
                    f.write(grace.format(system.output_name()))
            except:
                print("ion grace data failed")
                print(sys.exc_info())

            try:
                bulk_solv = self.bulk_solv_files[system].data
                energiesStr = ""
                for z, e in zip(coordinates, bulk_solv):
                    energiesStr += "{:f} {:f}\n".format(z, e)

                grace = graceScript.format(start_x, end_x, end_y, "bulk_solv", "{:s}") + energiesStr + "&"

                with open("grace/bulk_solv_" + system.output_name() + ".xgr", "w") as f:
                    f.write(grace.format(system.output_name()))
            except:
                print("bulk_solv grace data failed")
                print(sys.exc_info())

            try:
                bulk_ion = self.bulk_ion_files[system].data
                energiesStr = ""
                for z, e in zip(coordinates, bulk_ion):
                    energiesStr += "{:f} {:f}\n".format(z, e)

                grace = graceScript.format(start_x, end_x, end_y, "bulk_ion", "{:s}") + energiesStr + "&"

                with open("grace/bulk_ion_" + system.output_name() + ".xgr", "w") as f:
                    f.write(grace.format(system.output_name()))
            except:
                print("bulk_ion grace data failed")
                print(sys.exc_info())

    def paired_stats(self, sys1, sys2):
        diff_pmf = []
        diff_rdf = []
        diff_ion = []
        for a, d in zip(sys1, sys2):
            d_pmf = self.pmf_files[a].max - self.pmf_files[d].max
            diff_pmf.append(abs(d_pmf))

            d_rdf = self.solv_files[a].min - self.solv_files[d].min
            diff_rdf.append(abs(d_rdf))

            d_ion = self.ion_files[a].min - self.ion_files[d].min
            diff_ion.append(abs(d_ion))

        n = len(diff_pmf)

        if n == 0:
            return 0, ()

        mean_pmf = sum(diff_pmf) / n # average difference in pmf
        var_pmf = sum([(x - mean_pmf)**2 for x in diff_pmf]) / (n - 1)
        t_pmf = mean_pmf / (sqrt(var_pmf) / sqrt(n))

        mean_rdf = sum(diff_rdf) / n
        var_rdf = sum([(x - mean_rdf)**2 for x in diff_rdf]) / (n - 1)
        t_rdf = mean_rdf / (sqrt(var_rdf) / sqrt(n))

        mean_ion = sum(diff_ion) / n
        var_ion = sum([(x - mean_ion)**2 for x in diff_ion]) / (n - 1)
        t_ion = mean_ion / (sqrt(var_ion) / sqrt(n))

        return n, (mean_pmf, mean_rdf, mean_ion), (t_pmf, t_rdf, t_ion)


    def run(self, inp, grace = False, dipole = False):

        
        n_windows = 60
        idx = pd.IndexSlice

        data_types = ['solv', 'ion']
        col_indices = [('bulk', i) for i in data_types]
        col_indices.extend([(window, i) for window in range(n_windows) for i in ['z', 'pmf'] + data_types])

        row_index = pd.MultiIndex.from_product([[1,2], [7,10,14], ['acn','dce','h2o'], ['bmim','tma_tmea'], ['BF4', 'TMA', 'TMEA']])
        col_index = pd.MultiIndex.from_tuples(col_indices)

        s = pd.DataFrame(index=row_index, columns=col_index)

        print(s)
        print(s.dtypes)

        for key, data in self.pmf_files.items():
            pd_row = key.pandas_tuple()
            if len(data.data) != 60:
                continue
            s.loc[idx[pd_row], idx[:, 'z']] = data.data.keys()
            s.loc[idx[pd_row], idx[:, 'pmf']] = data.data.values()

        for key, data in self.solv_files.items():
            pd_row = key.pandas_tuple()
            if len(data.data) != 60:
                continue
            s.loc[idx[pd_row], idx[list(range(n_windows)), 'solv']] = data.data

        for key, data in self.ion_files.items():
            pd_row = key.pandas_tuple()
            if len(data.data) != 60:
                continue
            s.loc[idx[pd_row], idx[list(range(n_windows)), 'ion']] = data.data

        for key, data in self.bulk_solv_files.items():
            pd_row = key.pandas_tuple()
            s.loc[idx[pd_row], idx['bulk', 'solv']] = data.average

        for key, data in self.bulk_ion_files.items():
            pd_row = key.pandas_tuple()
            s.loc[idx[pd_row], idx['bulk', 'ion']] = data.average


        print(s)
        s.to_csv('data.csv')
        s.to_hdf('data.hd5', 'porous')






        return

        all_data = set(self.pmf_files.keys()) | set(self.solv_files.keys()) | set(self.ion_files.keys()) | set(self.bulk_solv_files.keys()) | set(self.bulk_ion_files.keys())
        queried_systems = find_sys(inp, all_data)

        for i in sorted(queried_systems):
            try:
                print("{:s}: pmf {:0.1f} {:0.1f} {:0.1f} solv {:0.1f} {:0.1f} {:0.1f} ion {:0.1f} {:0.1f} {:0.1f} bulk_solv {:0.1f} bulk_ion {:0.1f}".format(i,
                    self.pmf_files[i].start, self.pmf_files[i].max, self.pmf_files[i].end,
                    self.solv_files[i].start, self.solv_files[i].min, self.solv_files[i].end,
                    self.ion_files[i].start, self.ion_files[i].min, self.ion_files[i].end,
                    self.bulk_solv_files[i].average, self.bulk_ion_files[i].average
                    ))
            except:
                pass

        if grace:
            self.make_grace(all_data)
            return

        if dipole:
            for system in queried_systems:
                sys_path = system.sim_path()
                traj = sys_path + "/md_nvt_prod.dcd"
                top  = sys_path + "/md_nvt_prod.pdb"

                ion_pdb = "setup_data/" + (system.diffusingIon.upper() + ".pdb")

                bond_defs  = "ffdir/sapt_residues.xml"
                forcefield = "ffdir/sapt.xml"

                calc_dipole(system.diffusingIon, top, traj, bond_defs, forcefield, bulk = False, mol_pdb = ion_pdb)
            return

        ######
        #
        # run stats
        #
        ######

        systems_with_full_data = [x for x in queried_systems if x in self.pmf_files.keys() and x in self.solv_files.keys() and x in self.ion_files.keys()]

        for i in queried_systems:
            if i not in systems_with_full_data:
                print("no full data:", i)

        acn = []
        dce = []
        for i in systems_with_full_data:
            if i.solvent == "acn":
                analogue = None
                for j in queried_systems:
                    if j.solvent == "dce" and i.is_analogue(j):
                        analogue = j
                        break

                if analogue is not None:
                    acn.append(i)
                    dce.append(analogue)

        acn.sort()
        dce.sort()

        stats = self.paired_stats(acn, dce)
        if stats[0] > 1:
            print(stats[0], "pairs", stats[1], stats[2])
        else:
            print("no stats")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--filter", default="", help="filter systems")
    parser.add_argument("--write_grace", action="store_true", help="write grace data into a folder named grace")
    parser.add_argument("--dipole", action="store_true", help="print dipole information")
    args = parser.parse_args()

    dd = DiffusionData(sys.stdin.readlines())

    #
    # if being used as a single command
    #
    if len(sys.argv) > 1:
        dd.run(args.filter, args.write_grace, args.dipole)
        sys.exit(0)

    #
    # if being used as a REPL
    #
    sys.stdin = open("/dev/tty") # only works on Linux probably
    while True:
        inp = input()

        if inp == "q" or inp == "quit":
            break

