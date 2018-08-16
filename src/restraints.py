from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

import math
import numpy as np
from numpy import linalg, mod
from numpy.linalg import norm

from simulation import ionAtomNameMap

def findCenterOfPore(sim):
    """
    Find relevant C atoms on the edge of the pore.

    The center of the pore is defined by the center of the first and last
    sheets. Since the pore is defined across the boundary of the periodic
    image, and these two sheets are on the boundary, they constitute the
    center of the pore.

    The edges of the pore will be defined by the center of the two interior
    sheets. Since the umbrellas are moving in the positive z direction,
    the interior edge will be that defined by resid `len(sim.sheets) // 2`,
    and the exterior edge will by that defined by resid
    `len(sim.sheets) // 2 + 1`. Interior in this context means the sheets
    that face the solvent.

    It is important that the periodic image that is being used for the
    calculation contains the entire pore on one sheet. Using half of a sheet
    from consecutive periodic images will  will place the center of the
    pore at the boundary between sheets.

    Returns
    ______
    center of pore: [x0, y0, z0]
    interior edge of pore: [x0, y0, z0]
    exterior edge of pore: [x0, y0, z0]
    """

    # indices:
    # center: 0
    # interior: 1
    # exterior: 2
    coords = [np.array([0.0, 0.0, 0.0]) for x in range(3)]
    atomIndices = [[], [], []]

    graphRes = list(filter(lambda x: x.name == "grph", sim.simmd.topology.residues()))
    totalSheets = len(graphRes)

    # define the residue indices of the relevant sheets
    center1 = totalSheets // 2 - 1
    center2 = totalSheets // 2
    interior = totalSheets - 1
    exterior = 0

    # find pore center

    centerSheets = [graphRes[center1], graphRes[center2], graphRes[interior], graphRes[exterior]]
    # make sure that the resids we want match up to the objects we have
    assert centerSheets[0].index == center1
    assert centerSheets[1].index == center2
    assert centerSheets[2].index == interior
    assert centerSheets[3].index == exterior

    for i, sheet in enumerate(centerSheets):
        for atom in sheet._atoms:
            if atom.name in ["C354", "C274", "C148", "C68"]:
                if i <= 1:
                    atomIndices[0].append(atom.index)
                else:
                    atomIndices[i - 1].append(atom.index)

    # calculate centers
    state = sim.simmd.context.getState(getPositions=True, enforcePeriodicBox=True)
    prePotentialPositions = state.getPositions(asNumpy=True)
    boxDims = sim.simmd.topology.getPeriodicBoxVectors()
    box = np.mat(boxDims / nanometer)
    box_inv = linalg.inv(box)

    for i, center in enumerate(atomIndices):
        centerCoordList = []

        for index in center:
            nextAtomPosition = prePotentialPositions[index] / nanometer

            if len(centerCoordList) == 0:
                # this is the first atom in the group
                centerCoordList.append(nextAtomPosition)
                continue

            # make sure we are using the minimum image distance to the next atom from the first atom recorded

            dr = centerCoordList[0] - nextAtomPosition
            minimum_image_vector = minimum_image_distance(box, dr, box_inv)
            shift = minimum_image_vector - dr # really it is the negative of this, but we're just looking for non-zero terms so that doesn't matter

            if shift.any(): # if the shift is nonzero, a shift has occured
                print("Adjusting minimum distance for atom", index, "with shift",  shift)

            centerCoordList.append(centerCoordList[0] - minimum_image_vector)

        coords[i] = np.mean(centerCoordList, axis=0) # produce a vector of the averages wrt each coordinate

    return coords[0], coords[1], coords[2]


def minimum_image_distance(box, dr, box_inv = None):
    if box_inv is None:
        box_inv = linalg.inv(box)

    dr_box = box_inv.dot(dr).getA1() # getA1 flattens [[x,y,z]] (1x3 matrix) to [x,y,z] (vector)
    shift_box = np.array(list(map(lambda x: math.floor(x + 0.5), dr_box)))

    return dr - box.dot(shift_box).getA1()


def addIonUmbrellaPotential(sim, distanceFromPore, numbrella, kxy, ion_name):
    """
    Apply an umbrella potential window at the specified distanceFromPore
    from the edge of the pore for the Boron atom in the first BF4 residue.


    If the ion is on the wrong side of the system, it tries to force its
    way through the graphene sheet because the shortest path to the potential
    is through the sheet rather than going through the solvent.

    This is better than choosing the closest ion in the case that the closest
    ion of interest is on the wrong side of the sheet (this could occur with
    a dilute solvent, or other reasons).

    In order to avoid this, always start the ion at the center of the system.
    That way the shortest path to the potential will not cross through a graphene
    sheet.

    Parameters
    _________
    distanceFromPore: float
        Distance from the edge of the pore for the potential to exist.
    numbrella: int
        Number of umbrella potentials to be used
    kxy: float
        Force constant for the potential in the x, y dimensions.
    ion_name: string
        What residue pull through the pore

    Returns
    ______

    potentialCoordinates, [x0,y0,z0]: [float]
        Starting location of the potential
    dz: float
        Distance to move the umbrella after each iteration.
    """

    boxVecs = sim.simmd.context.getState().getPeriodicBoxVectors()
    boxCenter = np.array([0.0,0.0,0.0])
    for i in boxVecs:
        boxCenter += i / nanometer
    boxCenter /= 2

    print('Setting umbrella forces about equilibrated pore')
    centerCoords, interiorCoords, exteriorCoords = findCenterOfPore(sim)
    x0 = centerCoords[0]
    y0 = centerCoords[1]
    z0 = interiorCoords[2] - distanceFromPore
    print(10 * exteriorCoords, 10 * interiorCoords)

    rPore = exteriorCoords - interiorCoords
    periodicPoreLength = minimum_image_distance(np.mat(boxVecs / nanometer), rPore)
    if (rPore - periodicPoreLength).any():
        print("modified pore length from", rPore, "to", periodicPoreLength)
    totalDistanceToTravel = 2*distanceFromPore + abs(periodicPoreLength[2])
    dz = totalDistanceToTravel / numbrella

    # find an ion that is close to the pore, and is on the right side of the sheet
    aname = ionAtomNameMap[ion_name]
    ions_to_try = []
    for res in sim.simmd.topology.residues():
        if res.name == ion_name:
            for atom in res._atoms:
                if atom.name == aname:
                    ions_to_try.append((atom.index, res.id))

    index = None # just as a safeguard to make sure an error is always raised if this default index is used
    state = sim.simmd.context.getState(getPositions=True, enforcePeriodicBox=True)
    positions = state.getPositions()
    failure_iterations = 20
    for i in range(failure_iterations): # just in case there is not an ion near the starting location
                                        # we can iterate a few times to wait until ions have diffused
                                        # seems unlikely for this to be necessary
        for ion in ions_to_try:
            z = positions[ion[0]][2] / nanometer
            if z < interiorCoords[2] and z > boxCenter[2]:
                index = ion[0]
                resid = ion[1]
                break

        if index is not None:
            break
        else: # run the sim a little bit to move the ions around
            sim.simmd.step(100)

    print("Ion index:", index, "Resid:", resid)
    if index is None:
        raise IndexError("ion index was not found in the system")

    # Apply potential
    # PMF apply umbrella potential to z direction
    kz=2000.0  # kJ/mol/nm^2
    ZForce = CustomExternalForce("0.5*kz*periodicdistance(x,y,z,x,y,z0)^2")
    sim.system.addForce(ZForce)
    ZForce.addParticle(index)
    ZForce.addGlobalParameter('z0', z0)
    ZForce.addGlobalParameter('kz', kz)

    # PMF apply umbrella potential to x, y directions
    XYForce = CustomExternalForce("0.5*kxy*periodicdistance(x,y,z,x0,y0,z)^2")
    sim.system.addForce(XYForce)
    XYForce.addParticle(index)
    XYForce.addGlobalParameter('x0', x0)
    XYForce.addGlobalParameter('y0', y0)
    XYForce.addGlobalParameter('kxy', kxy)

    return [x0, y0, z0], dz, index

def restrainGrapheneSheets(sim, restrained_atoms, z0graph):
    state = sim.simmd.context.getState(getPositions=True)
    preGrapheneRestraintPositions = state.getPositions()

    all_atoms = [item for sublist in restrained_atoms for item in sublist]
    all_atoms = map(lambda x: x.index, all_atoms)
    for i, carbon in enumerate(all_atoms):
        x0graph = preGrapheneRestraintPositions[carbon][0] / nanometer
        y0graph = preGrapheneRestraintPositions[carbon][1] / nanometer

        # add restraining force to 10 atoms in the sheet
        kgraph=2000.0  # kJ/mol/nm^2
        GraphForce = CustomExternalForce("0.5*kgraph*periodicdistance(x,y,z, x0graph" + str(i) + ", y0graph" + str(i) + ", z0graph)^2")
        sim.system.addForce(GraphForce)
        GraphForce.addParticle(carbon)
        GraphForce.addGlobalParameter('kgraph', kgraph)
        GraphForce.addGlobalParameter('x0graph' + str(i), x0graph)
        GraphForce.addGlobalParameter('y0graph' + str(i), y0graph)
        GraphForce.addGlobalParameter('z0graph', z0graph)


def freezeSheets(sim, reset_mass = 0):
    graph_residues = list(filter(lambda r: r.name == "grph", sim.simmd.topology.residues()))
    carbon_mass = None
    for i, res in enumerate(graph_residues):
        for atom in res._atoms:
            if carbon_mass is None:
                carbon_mass = sim.system.getParticleMass(atom.index) / dalton

            sim.system.setParticleMass(atom.index, reset_mass * dalton)

    return carbon_mass

def bondGrapheneSheets(sim, numSheets):
    graph_residues = list(filter(lambda r: r.name == "grph", sim.simmd.topology.residues()))
    print("Restrained atoms ({:d} residues):".format(len(list(graph_residues))))
    restrained_atoms = [[] for y in range(2 * numSheets)]

    for i, res in enumerate(graph_residues):
        for atom in res._atoms:
            if sim.system.getParticleMass(atom.index) / dalton == 12.0108: # this mass was explicitly changed in the forcefield for atoms on the pore; it's just a marker for the pore-atoms
                restrained_atoms[i].append(atom)

    i = 0
    print(len(restrained_atoms[1]))
    while i < len(restrained_atoms):
        for j in range(len(restrained_atoms[i])):
            # make a bond between this sheet and the adjacent one
            print(restrained_atoms[i][j], restrained_atoms[i + 1][j])

            new_bond = sim.harmonicBondForce.addBond(restrained_atoms[i][j].index, restrained_atoms[i + 1][j].index, length=3.4 * angstrom, k=459403.2)

        # skip to the next pair of adjacent sheets
        i += 2

    if len(restrained_atoms) > 2:
        for i in range(len(restrained_atoms[0])):
            # make a bond between the first sheet and the last one
            sim.simmd.topology.addBond(restrained_atoms[0][i], restrained_atoms[-1][i])

    return restrained_atoms
