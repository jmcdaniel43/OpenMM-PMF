#!/usr/bin/env python

from __future__ import print_function
import sys
import re
import fileinput
from subprocess import call

A5  = [362, 361, 351, 261, 251, 250]
A7  = A5 +  [372, 352, 350, 262, 260, 240]
A10 = A7 +  [360, 259, 249, 340, 341, 241, 252, 363, 373, 272, 271, 371]
A14 = A10 + [38, 48, 58, 139, 149, 159, 169, 230, 231, 239, 242, 253, 263, 269, 270, 273, 281, 282, 283, 330, 331, 342, 342, 353, 364, 370, 374, 381, 382, 383, 384]



def makePore(pore, filenames):
    if pore == "5":
        atomlist = A5
    elif pore == "7":
        atomlist = A7
    elif pore == "10":
        atomlist = A10
    elif pore == "14":
        atomlist = A14
    else:
        raise ValueError("pore not one of <5, 7, 10, 14>")

    atomlist.sort()

    for file in filenames:
        for line in open(file).readlines():
            cont = False
            extension = file.split('.')[-1]

            for atom in atomlist:
                if re.search(r"C{0:d}[ \"]".format(atom), line, flags=re.UNICODE) is not None:
                    cont = True
                    break

            if cont:
                continue
            else:
                print(line, end="")
        print() # between files there is a newline

def removeDrudes(filenames):
    for f in filenames:
        call(["sed", "-i", "s/grp /grph/g", f])
        call(["sed", "-i", "s/TME /TMEA/g", f])
        call(["sed", "-i", "s/acn /acnt/g", f])
        call(["sed", "-i", "s/dch /dchl/g", f])
        for line in open(f).readlines():
            if re.search(' D[A-Z]', line) is not None:
                continue
            elif re.search('CONECT', line) is not None:
                continue
            else:
                print(line, end="")

if __name__ == "__main__":
    if len(sys.argv) > 2:
        if sys.argv[1] == "pore":
            makePore(sys.argv[2], sys.argv[3:])
        if sys.argv[1] == "drudes":
            removeDrudes(sys.argv[2:])
    else:
        print("USAGE: " + sys.argv[0] + " <pore [pore size: 5, 7, 10], drudes> <file names>")
