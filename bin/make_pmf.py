#!/usr/bin/env python

from multiprocessing import Pool
import fileinput
import os
import re
from os import path
from subprocess import call
from convert_umbrella_output import convert_umbrella_output
from simpath2outputname import simpath2outputname

whaminputFile = "whaminput"
whaminputRegex = re.compile("umbrella_([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+)\s+([0-9]+\.[0-9]+)")

startingDir = os.getcwd()

# for sys in map(lambda x: x.strip(), fileinput.input()):

def make_pmf(sys):
    os.chdir(startingDir)

    if not path.exists(sys):
        print("no system exists:", sys)
        return

    os.chdir(sys)

    try:
        os.mkdir("pmf")
    except FileExistsError:
        if path.exists("pmf/pmf"):
            print("pmf exists:", sys)
            # continue

    os.chdir("pmf")

    log_name = startingDir + "/output_logs/" + simpath2outputname([sys])[0] + ".log"

    try:
        whaminput = convert_umbrella_output(open(log_name), 2000, 10000, 60)
    except:
        print("error converting log:", log_name)
        return

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

    # wham hist_min hist_max num_bins tol temperature numpad metadatafile freefile num_MC_trials randSeed

    call(['wham', startWindow, endWindow, "60", "0.01", "300.0", "0", whaminputFile, "pmf", "1000", "143289"])

with Pool(12) as p:
    p.map(make_pmf, map(lambda x: x.strip(), fileinput.input()))
