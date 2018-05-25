#!/usr/bin/env python

import sys
import fileinput

from common import DiffusionSystem

def outputname2simpath(files):
    return_strings = []
    
    for line in fileinput.input():
        try:
            sys = DiffusionSystem(line)
        except ValueError as e:
            print(e)
            continue
        
        return_strings.append(sys.sim_path())

    return return_strings

if __name__ == '__main__':
    for i in outputname2simpath(fileinput.input()):
        print(i)
