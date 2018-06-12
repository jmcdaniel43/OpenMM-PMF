pmfFiles=`find $1 | grep --color=none ".*diff/ion_coordination.dat$"`

declare -a pmfData=()

for i in $pmfFiles
do
    nums=`cat $i | python -c "
from __future__ import print_function
import fileinput

nums = []
for i in fileinput.input():
    nums.append(float(i.strip()))
print(nums[0], min(nums), nums[-1])
"`
    pmfData=("${pmfData[@]}" "`echo $i | bin/simpath2outputname.py` $nums")

done

echo ${pmfData[@]} | python -c "
from __future__ import print_function
import fileinput
import sys

pmfMap = {}


iter = iter(fileinput.input()[0].split(' '))

pmfData = zip(iter, iter, iter, iter) # this groups the list into 4's; the iterator is mutably incremented each access
for i in pmfData:
    pmfMap[i[0]] = list(map(float, i[1:]))

for key in sorted(pmfMap.keys()):
    print(key, pmfMap[key][1], pmfMap[key][1], pmfMap[key][2])
"
