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

pmf_regex = re.compile("(\d+\.\d+)\s+(\d+\.\d+)\s+.*")

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


    def run(self, inp):
        args = inp.split(' ')
        attr = None
        if len(args) > 1:
            attr = args[1]

        systems_with_full_data = [x for x in self.pmf_files.keys() if x in self.rdf_files.keys() and x in self.ion_files.keys()]

        for i in (list(self.pmf_files.keys()) + list(self.rdf_files.keys())):
            if i not in systems_with_full_data:
                print("not included:", i)

        queried_systems = find_sys(inp, systems_with_full_data)
        for i in sorted(queried_systems):
            if attr is not None:
                # print("{:s} {:s} {:f} {:f} {:f}".format(i, attr)
                pass

            print("{:s}: pmf {:f} {:f} {:f} solv {:f} {:f} {:f} ion {:f} {:f} {:f} bulk_solv {:f} bulk_ion {:f}".format(i,
                self.pmf_files[i].start, self.pmf_files[i].max, self.pmf_files[i].end,
                self.rdf_files[i].start, self.rdf_files[i].min, self.rdf_files[i].end,
                self.ion_files[i].start, self.ion_files[i].min, self.ion_files[i].end,
                self.bulk_solv_files[i].average, self.bulk_ion_files[i].average
                ))

        acn = []
        dce = []
        for i in queried_systems:
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

if __name__ == '__main__':
    dd = DiffusionData(sys.stdin.readlines())

    #
    # if being used as a single command
    #
    if len(sys.argv) > 1:
        dd.run(sys.argv[1])
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
