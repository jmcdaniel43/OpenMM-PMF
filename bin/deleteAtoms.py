#!/usr/bin/env python

from __future__ import print_function
import sys
import re
import fileinput
from subprocess import call

A5  = list(map(str, [362, 361, 351, 261, 251, 250]))
A7  = list(map(str, [372, 352, 350, 262, 260, 240]))
A10 = list(map(str, [360, 259, 249, 340, 341, 241, 252, 363, 373, 272, 271, 371]))
A14 = list(map(str, [38, 48, 58, 139, 149, 159, 169, 230, 231, 239, 242, 253, 263, 269, 270, 273, 281, 282, 283, 330, 331, 342, 342, 353, 364, 370, 374, 381, 382, 383, 384]))
A17 = list(map(str, [27, 28, 37, 47, 57, 68, 128, 129, 138, 148, 158, 168, 179, 220, 221, 229, 232, 243, 254, 264, 274, 279, 280, 284, 320, 321, 332, 343, 354, 365, 375, 380, 385, 691, 692, 693, 694, 791, 792, 793, 794, 795]))

atomlists = [A5, A7, A10, A14, A17]



def makePore(pore, filenames):
    if pore == "5":
        a_idx = 0
    elif pore == "7":
        a_idx = 1
    elif pore == "10":
        a_idx = 2
    elif pore == "14":
        a_idx = 3
    else:
        raise ValueError("pore not one of <5, 7, 10, 14>")

    all_pore_atoms = [atom for sublist in atomlists[:a_idx + 1] for atom in sublist]

    atom_regex = re.compile("C(\d+)")

    for line in sys.stdin.readlines():
        match = atom_regex.findall(line)
        if len(match) == 0:
            print(line, end="")

        else:
            pore_atom = False
            for i in match:
                if i in all_pore_atoms: # if the atom is in the remove-list for each pore, dont re-print
                    pore_atom = True
                    break

            if pore_atom:
                continue

            # the above lines take care of the Bond attributes where multiple atoms are listed on a single line
            # from here on, only a single match is worth looking at because we are dealing with forcefield parameters
            elif match[0] in atomlists[a_idx + 1]:
                # if the atom is on the edge of the pore after removal (i.e., on the remove-list for the next pore)
                # change the atom type in the forcefield
                # in the PDB file, this command will have no effect (Cgt doesn't show up in the PDB)
                reassigned_atom_line = line.replace("Cgt", "Cpore").replace("12.011", "12.0108")
                print(reassigned_atom_line, end="")

            else:
                print(line, end="")

def removeDrudes(filenames):
    for f in filenames:
        call(["sed", "-i", "s/grp /grph/g", f])
        call(["sed", "-i", "s/TME /TMEA/g", f])
        call(["sed", "-i", "s/acn /acnt/g", f])
        call(["sed", "-i", "s/dch /dchl/g", f])
        call(["sed", "-i", "s/BMI /BMIM/g", f])
        for line in open(f).readlines():
            if re.search(' D[A-Z]', line) is not None:
                continue
            if re.search(' [A-Z]D', line) is not None:
                continue
            if re.search(' M ', line) is not None:
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
        print("USAGE: " + sys.argv[0] + " <pore [pore size: 5, 7, 10, 14], drudes> <file names>")
