import numpy as np
import math

def minimum_image_disp(box, dr, box_inv = None):
    """
    Computes the minimum image displacement for the given vector, dr.
    """
    if box_inv is None:
        box_inv = linalg.inv(box)

    dr_box = box_inv.dot(dr).getA1() # getA1 flattens [[x,y,z]] (1x3 matrix) to [x,y,z] (vector)
    shift_box = map(lambda x: math.floor(x + 0.5), dr_box)

    return dr - box.dot(list(shift_box)).getA1()

