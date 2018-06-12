import re
import sys
from pathlib import Path
from math import sqrt
from multiprocessing import Pool

from common import DiffusionSystem
from make_pmf import make_pmf
# from make_rdf import make_rdf

"""
Retrieves data from PMF and RDF files for statistical calculations
"""

pmf_files = {}
rdf_files = {}

pmf_regex = re.compile("(\d+\.\d+)\s+(\d+\.\d+)\s+.*")

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

    def __len__(self):
        return len(self.data)

def get_pmf(datapath):
    pmf_path = Path(datapath) / "pmf" / "pmf"
    if not pmf_path.exists():
        print("no such pmf files", datapath, file=sys.stderr)
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

        for k, v in pmf_files.items():
            if len(v) != len(pmf_values):
                print("this pmf file is not like the others", datapath, "(compared to", k, ")", file=sys.stderr)
            break

    if len(pmf_values) < 1:
        return

    system = DiffusionSystem(datapath)
    pmf_files[system.output_name()] = PMF_Data(pmf_values)

def get_rdf(datapath):
    rdf_path = Path(datapath) / "rdf.log"
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

        for k, v in rdf_files.items():
            if len(v) != len(rdf_values):
                print("this rdf file is not like the others", datapath, "(compared to", k, ")", file=sys.stderr)
            break

    if len(rdf_values) < 1:
        return

    system = DiffusionSystem(datapath)
    rdf_files[system.output_name()] = RDF_Data(rdf_values)

def find_sys(inp, systems):
    try:
        system_filter = re.compile(inp)
    except:
        # probably an invalid regex
        print("regex excceptiong", file=sys.stderr)

    queried_systems = [x for x in systems if system_filter.match(x) is not None] 

    return_systems = []
    for i in queried_systems:
        return_systems.append("{:s}: pmf {:f} {:f} {:f} rdf {:f} {:f} {:f}".format(i,
            pmf_files[i].start, pmf_files[i].max, pmf_files[i].end,
            rdf_files[i].start, rdf_files[i].min, rdf_files[i].end,
            ))

    return_systems.sort()

    if len(return_systems) < 1:
        return []

    return return_systems

def paired_stats(systems):
    acn = []
    dce = []
    for i in map(DiffusionSystem, systems):
        if i.solvent == "acn":
            acn.append(i)
        else:
            dce.append(i)

    acn.sort()
    dce.sort()

    diff_pmf = []
    diff_rdf = []
    for a, d in zip(acn, dce):
        d_pmf = pmf_files[a.output_name()].max - pmf_files[d.output_name()].max
        diff_pmf.append(abs(d_pmf))

        d_rdf = rdf_files[a.output_name()].min - rdf_files[d.output_name()].min
        diff_rdf.append(abs(d_rdf))

    n = len(diff_rdf)

    mean_pmf = sum(diff_pmf) / n
    stdev_pmf = sum([(x - mean_pmf)**2 for x in diff_pmf]) / (n - 1)
    t_pmf = mean_pmf / (stdev_pmf / sqrt(n))

    mean_rdf = sum(diff_rdf) / n
    stdev_rdf = sum([(x - mean_rdf)**2 for x in diff_rdf]) / (n - 1)
    t_rdf = mean_rdf / (stdev_rdf / sqrt(n))

    return n, (t_pmf, t_rdf)


def main(inp):
    systems_with_full_data = [x for x in pmf_files.keys() if x in rdf_files.keys()]

    queried_systems = find_sys(inp, systems_with_full_data)
    for i in queried_systems:
        print(i)

    stats = paired_stats(queried_systems)
    print(stats[0], stats[1])

#     first = DiffusionSystem(queried_systems[0])
#     same_system = True
#     for i in queried_systems:
#         j = DiffusionSystem(i)
#         if j.tmaPresent != first.tmaPresent:
#             same_system = False
# 
#     if same_system:
#         if not first.tmaPresent: # add newline spacers to fit the excel sheet
#             l = ["" for x in range(2 * len(return_systems))]
#             l[::2] = return_systems
#             return_systems = l


def get_data(path):
    get_pmf(path)
    get_rdf(path)


if __name__ == '__main__':
    for line in sys.stdin.readlines():
        path = line.strip()
        get_pmf(path)
        get_rdf(path)

    #
    # if being used as a single command
    #
    if len(sys.argv) > 1:
        main(sys.argv[1])
        sys.exit(0)

    #
    # if being used as a REPL
    #
    sys.stdin = open("/dev/tty") # only works on Linux probably
    while True:
        inp = input()

        if inp == "q" or inp == "quit":
            break

        main(inp)
