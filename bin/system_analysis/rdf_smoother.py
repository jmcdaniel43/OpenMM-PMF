#!/usr/bin/env python

import sys

smooth_range = 5
smooth_limit = smooth_range // 2

def smooth(rdf_data):
    assert len(rdf_data) >= smooth_range
    retstring = []

    for i in range(len(rdf_data)):
        if i - smooth_limit < 0:
            nearest_neighbors = rdf_data[:smooth_range]
        elif i + smooth_limit >= len(rdf_data):
            nearest_neighbors = rdf_data[len(rdf_data) - smooth_range:]
        else:
            nearest_neighbors = rdf_data[i - smooth_limit:i + smooth_limit + 1]

        nearest_neighbors = map(float, nearest_neighbors)
        smoothed = sum(nearest_neighbors) / smooth_range
        retstring.append(smoothed)

    return retstring

if __name__ == '__main__':
    for i in smooth(sys.stdin.readlines()):
        print(i)
