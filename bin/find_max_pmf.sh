pmfFiles=`find $1 | grep --color=none ".*diff/pmf/pmf$"`

declare -a pmfData=()

for i in $pmfFiles
do
    nums=`grep "^[0-9][0-9]\.[0-9]*" $i | awk '{print $2}' \
        | python -c "
import fileinput

nums = []
for i in fileinput.input():
    nums.append(float(i.strip()))
print(nums[0], max(nums), nums[-1])
"`
    pmfData=("${pmfData[@]}" "`echo $i | bin/simpath2outputname.py` $nums")

done

echo ${pmfData[@]} | python -c "
import fileinput
import sys

pmfMap = {}

iter = iter(fileinput.input()[0].split(' '))

pmfData = zip(iter, iter, iter, iter) # this groups the list into 4's; the iterator is mutably incremented each access
for i in pmfData:
    pmfMap[i[0]] = list(map(float, i[1:]))

for key in sorted(pmfMap.keys()):
    print(key, pmfMap[key][0], pmfMap[key][1], pmfMap[key][2])
"
