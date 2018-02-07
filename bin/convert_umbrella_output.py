#!/usr/bin/env python

import re
import sys
import glob

def main():
    if len(sys.argv) < 7:
        print("USAGE: <input file> <spring constant> <start> <dz> <nstep> <numbrella>")
        exit()

    springConstant = float(sys.argv[2])
    start          = float(sys.argv[3])
    dz             = float(sys.argv[4])
    nstep          = int(sys.argv[5])
    numbrella      = int(sys.argv[6])

    startRegex = re.compile("Starting Production")
    with open(sys.argv[1]) as f:
        lines = f.readlines()
        for line in range(len(lines)):
            if startRegex.match(lines[line]) is not None:
                line += 1
                break
        for i in range(numbrella):
            z = start + dz*(i + 1)
            with open("umbrella_{:0.15f}".format(z), "w") as ofile:
                for j in range(nstep):
                    pos  = lines[line].split(',')
                    zpos = float(pos[2].split(')')[0]) * 10 # convert nm to angstrom
                    ofile.write("{:d}   {:0.15f}\n".format(j, zpos))
                    line += 1

    for file in glob.glob("umbrella_*"):
        x0 = file[9:]
        k = springConstant / 4.184 / 100 # convert springConstant from kJ/mol/nm^2 to kcal/mol/A^2
        print("{:s}    {:s}      {:0.2f}".format(file, x0, k))
        


if __name__ == "__main__":
    main()
