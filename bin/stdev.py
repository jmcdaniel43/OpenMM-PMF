#!/usr/bin/python

import re
from glob import glob
from math import sqrt

z_regex = re.compile("\d+\s+(\d+\.\d+)")

def stdev(expected_mean, whaminput):
    zs = []
    for line in whaminput.splitlines():
        match = z_regex.match(line)
        if match is None:
            continue

        zs.append(float(match.group(1)))
            
    n = len(zs)
    return sqrt(sum([(x - expected_mean)**2 for x in zs]) / (n - 1))

for idx, i in enumerate(glob("umbrella*")):
    filename = i.split(" ")[0].strip()
    with open(filename) as f:
        expected_mean = float(filename.split("_")[1])
        cur = stdev(expected_mean, f.read())
        if cur > 1: # if the stdev is greater than 1 Angstrom away from the center of the potential
                    # then this window probably has something weird going on
            print("stdev is too high {:f} window {:s} (#{:d})".format(cur, filename, idx))
            # return
        else:
            print(cur)
