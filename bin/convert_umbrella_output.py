#!/usr/bin/env python

import re
import sys
import glob


zDimensionRegexString = r"(-?[0-9]+\.[0-9]+(e-[0-9]+)?),? (-?[0-9]+\.[0-9]+(e-[0-9]+)?),? (-?[0-9]+\.[0-9]+(e-[0-9]+)?)"
maxZDimensionRegex = re.compile("Box dimensions:  \(" + zDimensionRegexString + "\) nm")
startingPotentialDimensionRegex1 = re.compile("Center of umbrella potential \(nm\) \[" + zDimensionRegexString + "\]")
startingPotentialDimensionRegex2 = re.compile("Center of umbrella potential \(x,y,z \(nm\)\): " + zDimensionRegexString)
startingPotentialDimensionRegex3 = re.compile("Center of umbrella potential \(nm\): x0 y0 z0 " + zDimensionRegexString)
zLocationRegex = re.compile("\(" + zDimensionRegexString + "\)")
startRegex = re.compile("Starting Production")
dzRegex = re.compile("dz (\d+\.\d+)")

boundaryTolerance = 20

def convert_umbrella_output(fileHandle, springConstant, nstep, numbrella):
    lines = fileHandle.readlines()
    for line in range(len(lines)): # find start and maximum z dimension from output file
        maxZ_match = maxZDimensionRegex.match(lines[line])
        startZ_match1 = startingPotentialDimensionRegex1.match(lines[line])
        startZ_match2 = startingPotentialDimensionRegex2.match(lines[line])
        startZ_match3 = startingPotentialDimensionRegex3.match(lines[line])
        dz_match = dzRegex.match(lines[line])

        if maxZ_match is not None:
            maxZ = float(maxZ_match.group(5)) * 10 # convert nm to angstrom

        elif startZ_match1 is not None:
            start = float(startZ_match1.group(5)) * 10 # convert nm to angstrom
        elif startZ_match2 is not None:
            start = float(startZ_match2.group(5)) * 10 # convert nm to angstrom
        elif startZ_match3 is not None:
            start = float(startZ_match3.group(5)) * 10 # convert nm to angstrom

        elif dz_match is not None:
            dz = float(dz_match.group(1)) * 10 # convert nm to angstrom

        elif startRegex.match(lines[line]) is not None:
            line += 1
            break

    try:
        dz
    except NameError:
        dz = 0.4 # default value of 0.4A for older simulations

    umbrellas = []

    # copy z location data into binned files
    for i in range(numbrella):
        # if False:
        #     z = start + dz * (i + 1) # some early simulations preincremented dz, so use this line
        #     print("starting at dz * (i + 1)", file=sys.stderr)
        z = start + dz * i # center of current umbrella potential

        umbrellas.append("umbrella_{:0.15f}".format(z))

        j = 0
        with open(umbrellas[-1], "w") as ofile:
            while j < nstep:
                zpos_match = zLocationRegex.match(lines[line])
                if zpos_match is not None:
                    zpos = float(zpos_match.group(5)) * 10 # convert nm to angstrom
                    if abs(zpos - z) > boundaryTolerance:   # then we will assume the ion has crossed the boundary of the period image
                        zpos = maxZ + zpos                  # just add into +z dimension past the end of the box
                                                            # this prevents us from having to make wham deal with wrapping
                    ofile.write("{:d}   {:0.15f}\n".format(j, zpos))
                    j += 1

                # else: # this is most likely a marker for ion movement
                    # j = j : don't count this line
                    # print(lines[line], file=sys.stderr, end="")

                line += 1

    return_string = ""
    for file in umbrellas:
        x0 = file[9:]
        k = springConstant / 4.184 / 100 # convert springConstant from kJ/mol/nm^2 to kcal/mol/A^2
        return_string += "{:s}    {:s}      {:0.2f}\n".format(file, x0, k)

    return return_string
        
if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("USAGE: <input file> <spring constant> <nstep> <numbrella>")
        exit()

    springConstant = float(sys.argv[2])
    nstep          = int(sys.argv[3])
    numbrella      = int(sys.argv[4])

    with open(sys.argv[1]) as f:
        print(convert_umbrella_output(f, springConstant, nstep, numbrella))
