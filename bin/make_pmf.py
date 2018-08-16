#!/usr/bin/env python

"""
Must be run in the root directory of the project, so that this script has access
to the ouput logs from the simulations.

"""

from multiprocessing import Pool
import os
import sys
import re
import traceback
import argparse
from os import path
from subprocess import call
from convert_umbrella_output import convert_umbrella_output
from simpath2outputname import simpath2outputname
from math import sqrt

whaminputFile = "whaminput"
whaminputRegex = re.compile("umbrella_([0-9]+\.[0-9]+)\s+(-?[0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+)")
z_regex = re.compile("\d+\s+(-?\d+\.\d+)")

startingDir = os.getcwd()

def stdev(expected_mean, whaminput):
    zs = []
    for line in whaminput.splitlines():
        match = z_regex.match(line)
        if match is None:
            continue

        zs.append(float(match.group(1)))
            
    n = len(zs)
    return sqrt(sum([(x - expected_mean)**2 for x in zs]) / (n - 1))

def make_pmf(system, force=False):
    os.chdir(startingDir)

    ####
    #
    # logistics for finding pmf directory
    #
    ####

    if not path.exists(system):
        print("no system exists:", system)
        return False

    try:
        os.chdir(system)
        os.mkdir("pmf")
    except NotADirectoryError:
        print("system not present:", system)
        return False
    except FileExistsError:
        if path.exists("pmf/pmf"):
            print("pmf exists:", system, file=sys.stderr)
            if not force:
                return
            print("overwriting data")

    os.chdir("pmf")

    ####
    #
    # start calculating pmf
    #
    ####

    log_name = startingDir + "/output_logs/" + simpath2outputname(system) + ".log"

    try:
        whaminput, pore_center = convert_umbrella_output(open(log_name), 2000, 60)
    except:
        print("error converting log:", log_name)
        traceback.print_exc()
        return False

    # check to make sure none of the windows have wild deviations from the expected location of the potential
    for idx, i in enumerate(whaminput.splitlines()):
        filename = i.split(" ")[0]
        with open(filename) as f:
            expected_mean = float(filename.split("_")[1])
            cur = stdev(expected_mean - pore_center, f.read())
            if cur > 1: # if the stdev is greater than 1 Angstrom away from the center of the potential
                        # then this window probably has something weird going on
                print("{:s} stdev is too high {:f} window {:s} (#{:d})".format(system, cur, filename, idx))

    with open(whaminputFile, "w") as f:
        f.write(whaminput)

    # the first line of this should be the starting window
    # the last line will be the final window
    # since wham needs a starting and ending range for histogramming, we'll take these
    lines = whaminput.splitlines()

    wham1 = whaminputRegex.match(lines[0])
    if wham1 is not None:
        startWindow = wham1.group(2)

    wham2 = whaminputRegex.match(lines[-1])
    if wham2 is not None:
        endWindow = wham2.group(2)

    # the parameters in the following call have the descriptions below
    # wham hist_min hist_max num_bins tol temperature numpad metadatafile freefile num_MC_trials randSeed

    with open("whamout", "w") as pipeOut:
        call(['wham', startWindow, endWindow, "60", "0.01", "300.0", "0", whaminputFile, "pmf", "1000", "143289"], stdout=pipeOut)

    print("pmf calculated:", system)

    # make sure we end up in the directory where we started
    os.chdir(startingDir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--forceOverwrite", action="store_true", help="force overwrite of existing data")
    args = parser.parse_args()

    systems = []
    for line in sys.stdin.readlines():
        systems.append(line.strip())

    def force_make_pmf(system):
        make_pmf(system, force=True)

    with Pool(12) as p:
        if args.forceOverwrite:
            p.map(force_make_pmf, systems)
        else:
            p.map(make_pmf, systems)
