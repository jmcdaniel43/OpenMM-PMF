#!/usr/bin/env python

import sys
import fileinput

from common import DiffusionSystem

def simpath2outputname(line):
    sys = DiffusionSystem(line)
    return sys.output_name()

if __name__ == '__main__':
    for i in fileinput.input():
        try:
            print(simpath2outputname(i))
        except ValueError as e:
            print(e)
