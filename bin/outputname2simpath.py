#!/usr/bin/env python

import sys
import fileinput

from common import DiffusionSystem

def outputname2simpath(line):
    sys = DiffusionSystem(line)
    return sys.sim_path()

if __name__ == '__main__':
    for i in fileinput.input():
        try:
            print(outputname2simpath(i))
        except ValueError as e:
            print(e)
