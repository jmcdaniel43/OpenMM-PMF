#!/usr/bin/env python

import sys
from math import log
import matplotlib.pyplot as plt

nums = []
for line in sys.stdin.readlines():
    nums.append(float(line.strip()))

# plt.plot(nums)
plt.hist(nums, 12)
plt.show()
