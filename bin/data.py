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
        print("no such pmf file", pmf_path, file=sys.stderr)
        print("attempting to calculate", file=sys.stderr)
        if make_pmf(datapath) is not None:
            print("pmf calculation failed", file=sys.stderr)
            return

    pmf_values = {}

    with open(pmf_path) as pmf_file:
        for line in pmf_file:
            pmf_value = pmf_regex.match(line)
            if pmf_value is not None:
                pmf_values[float(pmf_value.group(1))] = float(pmf_value.group(2))

    if len(pmf_values) < 1:
        return

    return PMF_Data(pmf_values)

def get_rdf(datapath, ion = False, bulk = False):
    rdf_path = Path(datapath)
    if not rdf_path.exists():
        print("no such rdf file", rdf_path, file=sys.stderr)

        print("attempting to calculate", file=sys.stderr)
        if make_rdf(datapath, rdf_path, ion = ion, bulk = bulk) is not None:
            print("rdf calculation failed", file=sys.stderr)
            return

    rdf_values = []

    with open(rdf_path) as rdf_file:
        for line in rdf_file:
            rdf_values.append(float(line.strip()))

    if len(rdf_values) < 1:
        return

    return RDF_Data(rdf_values)

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

            for ion in [True, False]:
                for bulk in [True, False]:
                    if bulk:
                        if ion:
                            extension = "bulk_ion"
                            files = self.bulk_ion_files
                        else:
                            extension = "bulk_solvation"
                            files = self.bulk_solv_files
                    else:
                        if ion:
                            extension = "ion_coordination"
                            files = self.ion_files
                        else:
                            extension = "solv"
                            files = self.solv_files

                    rdf_path = Path(system.sim_path() + "/" + extension + ".dat")
                    if not rdf_path.exists():
                        print("no such rdf file", rdf_path, file=sys.stderr)

                        print("attempting to calculate", file=sys.stderr)
                        rdf_values = make_rdf(system.sim_path(), extension, ion = ion, bulk = bulk)
                    else:
                        rdf_values = []

                        with open(rdf_path) as rdf_file:
                            for line in rdf_file:
                                rdf_values.append(float(line.strip()))

                    if len(rdf_values) < 1:
                        print("error: rdf_values has length 0")
                        continue

                    rdf = RDF_Data(rdf_values)

                    if rdf is not None:
                        for k, v in files.items():
                            if len(v) != len(rdf):
                                print("this rdf file is not like the others", path, "(compared to", k, ")", file=sys.stderr)
                            break

                        files[system] = rdf

    def run(self, label = 'porous'):
        n_windows = 60
        idx = pd.IndexSlice

        data_types = ['solv', 'ion']
        col_indices = [('bulk', i) for i in data_types]
        col_indices.extend([('nonbulk', i) for i in ['z', 'pmf', 'solv', 'ion']])

        row_index = pd.MultiIndex.from_product([[1,2], [7,10,14], ['acn','dce'], ['bmim','tma_tmea'], ['BF4', 'TMA', 'TMEA']])
        col_index = pd.MultiIndex.from_tuples(col_indices)

        s = pd.DataFrame(index=row_index, columns=col_index)
        s = s.astype(object)
        print(s.dtypes)

        for key, data in self.pmf_files.items():
            if key.solvent == 'h2o':
                continue
            pd_row = key.pandas_tuple()

            if len(data.data) != 60:
                continue
            s.loc[idx[pd_row], idx['nonbulk', 'z']] = list(data.data.keys())
            s.loc[idx[pd_row], idx['nonbulk', 'pmf']] = list(data.data.values())

        for key, data in self.solv_files.items():
            if key.solvent == 'h2o':
                continue
            pd_row = key.pandas_tuple()
            if len(data.data) != 60:
                continue
            s.loc[idx[pd_row], idx['nonbulk', 'solv']] = list(data.data)

        for key, data in self.ion_files.items():
            if key.solvent == 'h2o':
                continue
            pd_row = key.pandas_tuple()
            if len(data.data) != 60:
                continue
            s.loc[idx[pd_row], idx['nonbulk', 'ion']] = list(data.data)

        for key, data in self.bulk_solv_files.items():
            if key.solvent == 'h2o':
                continue
            pd_row = key.pandas_tuple()
            s.loc[idx[pd_row], idx['bulk', 'solv']] = data.average

        for key, data in self.bulk_ion_files.items():
            if key.solvent == 'h2o':
                continue
            pd_row = key.pandas_tuple()
            s.loc[idx[pd_row], idx['bulk', 'ion']] = data.average


        print(s)
        s.to_csv('data_new.csv')
        s.to_hdf('data_new.hd5', label)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--label", default='porous')
    parser.add_argument("--force", action='store_true')
    args = parser.parse_args()

    dd = DiffusionData(sys.stdin.readlines())

    dd.run(args.label)
