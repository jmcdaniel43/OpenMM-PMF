#!/usr/bin/env python

import sys
import os
from subprocess import call

cwd = os.getcwd()

for sysPath in sys.argv[1:]:
    os.chdir(cwd)
    os.chdir(sysPath)
