import re
import sys
from pathlib import Path
from math import sqrt
from multiprocessing import Pool
import argparse

from common import DiffusionSystem
from make_pmf import make_pmf
from make_grace import graceScript
# from make_rdf import make_rdf

"""
Retrieves data from PMF and RDF files for statistical calculations
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
        self.rdf_files = {}
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

            rdf = get_rdf(path + "/rdf.dat")
            if rdf is not None:
                for k, v in self.rdf_files.items():
                    if len(v) != len(rdf):
                        print("this rdf file is not like the others", path, "(compared to", k, ")", file=sys.stderr)
                    break

                self.rdf_files[system] = rdf

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

    def paired_stats(self, sys1, sys2):
        diff_pmf = []
        diff_rdf = []
        diff_ion = []
        for a, d in zip(sys1, sys2):
            d_pmf = self.pmf_files[a].max - self.pmf_files[d].max
            diff_pmf.append(abs(d_pmf))

            d_rdf = self.rdf_files[a].min - self.rdf_files[d].min
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


    def run(self, inp, grace = False):
        all_data = set(self.pmf_files.keys()) | set(self.rdf_files.keys()) | set(self.ion_files.keys()) | set(self.bulk_solv_files.keys()) | set(self.bulk_ion_files.keys())
        queried_systems = find_sys(inp, all_data)
        for i in sorted(queried_systems):
            continue
            try:
                print("{:s}: pmf {:f} {:f} {:f} solv {:f} {:f} {:f} ion {:f} {:f} {:f} bulk_solv {:f} bulk_ion {:f}".format(i,
                    self.pmf_files[i].start, self.pmf_files[i].max, self.pmf_files[i].end,
                    self.rdf_files[i].start, self.rdf_files[i].min, self.rdf_files[i].end,
                    self.ion_files[i].start, self.ion_files[i].min, self.ion_files[i].end,
                    self.bulk_solv_files[i].average, self.bulk_ion_files[i].average
                    ))
            except:
                pass

        if grace:
            for system in all_data:
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
                    solv = self.rdf_files[system].data
                    energiesStr = ""
                    for z, e in zip(coordinates, solv):
                        energiesStr += "{:f} {:f}\n".format(z, e)

                    grace = graceScript.format(start_x, end_x, end_y, "solv", "{:s}") + energiesStr + "&"

                    with open("grace/solv_" + system.output_name() + ".xgr", "w") as f:
                        f.write(grace.format(system.output_name()))
                except:
                    print("rdf grace data failed")
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


        ######

        ######

        ######
        #
        # run stats
        #
        ######
        systems_with_full_data = [x for x in queried_systems if x in self.pmf_files.keys() and x in self.rdf_files.keys() and x in self.ion_files.keys()]
        for i in (list(self.pmf_files.keys()) + list(self.rdf_files.keys())):
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
    args = parser.parse_args()

    dd = DiffusionData(sys.stdin.readlines())

    #
    # if being used as a single command
    #
    if len(sys.argv) > 1:
        dd.run(args.filter, args.write_grace)
        sys.exit(0)

    #
    # if being used as a REPL
    #
    sys.stdin = open("/dev/tty") # only works on Linux probably
    while True:
        inp = input()

        if inp == "q" or inp == "quit":
            break

        dd.run(inp)
