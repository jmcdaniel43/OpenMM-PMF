#!/usr/bin/env python

import sys
import fileinput

from common import DiffusionSystem

def simpath2outputname(files):
    return_strings = []

    for line in files:
        try:
            sys = DiffusionSystem(line)
        except ValueError as e:
            print(e)
            continue
        
        return_strings.append(sys.output_name())

    if len(return_strings) == 0:
        raise ValueError("no simpaths could be converted")

    return return_strings

if __name__ == '__main__':
    for i in simpath2outputname(fileinput.input()):
        print(i)
