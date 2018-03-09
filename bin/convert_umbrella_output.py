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

    zDimensionRegexString = r"\((-?[0-9]+\.[0-9]+(e-[0-9]+)?), (-?[0-9]+\.[0-9]+(e-[0-9]+)?), (-?[0-9]+\.[0-9]+(e-[0-9]+)?)\) nm"
    maxZDimensionRegex = re.compile("Box dimensions:  " + zDimensionRegexString)
    zLocationRegex = re.compile(zDimensionRegexString)
    startRegex = re.compile("Starting Production")
    with open(sys.argv[1]) as f:
        # find start and maximum z dimension from output file
        lines = f.readlines()
        for line in range(len(lines)):
            z_match = maxZDimensionRegex.match(lines[line])
            if z_match is not None:
                maxZ = float(z_match.group(5)) * 10 # convert nm to angstrom
            if startRegex.match(lines[line]) is not None:
                line += 1
                break

        # copy z location data into binned files
        for i in range(numbrella):

            z = start + dz * i # center of current umbrella potential

            with open("umbrella_{:0.15f}".format(z), "w") as ofile:
                for j in range(nstep):
                    zpos_match = zLocationRegex.match(lines[line])
                    if zpos_match is not None:
                        zpos = float(zpos_match.group(5)) * 10 # convert nm to angstrom
                        if zpos - start < 0: # then the ion has crossed the boundary of the period image
                            zpos = maxZ + zpos # just add into +z dimension past the end of the box
                                # this prevents us from having to make wham deal with wrapping
                    else: # this is probably if the file is too short, wrong parameters
                        print(lines[line])

                    # pos  = lines[line].split(',')
                    # try:
                    #     zpos = float(pos[2].split(')')[0]) * 10 # convert nm to angstrom
                    #     if zpos - start < 0: # then the ion has crossed the boundary of the period image
                    #         zpos = maxZ + zpos # just add into +z dimension past the end of the box
                    #             # this prevents us from having to make wham deal with wrapping
                    #     print(pos)
                    ofile.write("{:d}   {:0.15f}\n".format(j, zpos))
                    line += 1

    for file in glob.glob("umbrella_*"):
        x0 = file[9:]
        k = springConstant / 4.184 / 100 # convert springConstant from kJ/mol/nm^2 to kcal/mol/A^2
        print("{:s}    {:s}      {:0.2f}".format(file, x0, k))
        


if __name__ == "__main__":
    main()
