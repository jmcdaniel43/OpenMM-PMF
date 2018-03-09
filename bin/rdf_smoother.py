#!/usr/bin/env python

import sys

smooth_range = 5
smooth_limit = smooth_range // 2

rdf_data = sys.stdin.readlines()
assert len(rdf_data) >= smooth_range

for i in range(len(rdf_data)):
    if i - smooth_limit < 0:
        nearest_neighbors = rdf_data[:smooth_range]
    elif i + smooth_limit >= len(rdf_data):
        nearest_neighbors = rdf_data[len(rdf_data) - smooth_range:]
    else:
        nearest_neighbors = rdf_data[i - smooth_limit:i + smooth_limit + 1]

    nearest_neighbors = map(float, nearest_neighbors)
    smoothed = sum(nearest_neighbors) / smooth_range
    print(smoothed)

