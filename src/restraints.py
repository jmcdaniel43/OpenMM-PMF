from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

import numpy as np

def findCenterOfPore(sim, poreSize):
    """
    Find relevant C atoms on the edge of the pore.

    The center of the pore is defined by the center of the first and last
    sheets. Since the pore is defined across the boundary of the periodic
    image, and these two sheets are on the boundary, they constitute the
    center of the pore.

    The edges of the pore will be defined by the center of the two interior
    sheets. Since the umbrellas are moving in the positive z direction,
    the interior edge will be that defined by resid len(sim.sheets) // 2,
    and the exterior edge will by that defined by resid
    len(sim.sheets) // 2 + 1.

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
    coords = [[0.0] * 3, [0.0] * 3, [0.0] * 3]
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

    centerSheets = [graphRes[center1], graphRes[center2],
            graphRes[interior], graphRes[exterior]]
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
    boxDims = sim.simmd.topology.getUnitCellDimensions()

    for i, center in enumerate(atomIndices):
        for index in center:
            nextAtomPosition = prePotentialPositions[index] / nanometer

            if coords[i].all() == 0.0: # not yet initialized
                coords[i] = nextAtomPosition / len(center)
            else:
                for j in range(len(coords[i])):
                    # we need to make sure that the other atom is on the same sheet image
                    if coords[i][j] - nextAtomPosition[j] > poreSize + 0.5:
                        # if the nextAtom is greater than 5 Angstroms more than the size of the pore
                        #  we need to subtract the size of the box in the current dimension
                        # to bring us to the correct atom position
                        print("Adjusting", j, "dimension for atom", index)
                        nextAtomPosition[j] -= boxDims[j] / nanometer

                coords[i] += nextAtomPosition / len(center)


    return coords[0], coords[1], coords[2]


def addIonUmbrellaPotential(sim, distanceFromPore, numbrella, kxy, index, poreSize):
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
    atomselection: int
        Index of the ion that the potential should operate on.

    Returns
    ______

    potentialCoordinates, [x0,y0,z0]: [float]
        Starting location of the potential
    dz: float
        Distance to move the umbrella after each iteration.
    """

    boxVecs = sim.simmd.context.getState().getPeriodicBoxVectors(asNumpy=True)
    boxCenter = np.array([0.0,0.0,0.0])
    for i in boxVecs:
        boxCenter += i / nanometer
    boxCenter /= 2
    print("Bringing ion to center of system:", boxCenter)

    print('Setting umbrella forces about equilibrated pore')
    centerCoords, interiorCoords, exteriorCoords = findCenterOfPore(sim, poreSize)
    x0 = centerCoords[0]
    y0 = centerCoords[1]
    z0 = interiorCoords[2] - distanceFromPore

    totalDistanceToTravel = distanceFromPore * 2 + abs(exteriorCoords[2] - interiorCoords[2])
    dz = totalDistanceToTravel / numbrella

    # Apply potential
    # PMF apply umbrella potential to z direction
    kz=2000.0  # kJ/mol/nm^2
    ZForce = CustomExternalForce("0.5*kz*periodicdistance(x,y,z,x,y,z0)^2")
    sim.system.addForce(ZForce)
    ZForce.addParticle(index)
    ZForce.addGlobalParameter('z0', boxCenter[2])
    ZForce.addGlobalParameter('kz', kz)

    # PMF apply umbrella potential to x, y directions
    XYForce = CustomExternalForce("0.5*kxy*periodicdistance(x,y,z,x0,y0,z)^2")
    sim.system.addForce(XYForce)
    XYForce.addParticle(index)
    XYForce.addGlobalParameter('x0', x0)
    XYForce.addGlobalParameter('y0', y0)
    XYForce.addGlobalParameter('kxy', kxy)

    state = sim.simmd.context.getState(getPositions=True, enforcePeriodicBox=True)
    sim.simmd.context.reinitialize()
    sim.simmd.context.setPositions(state.getPositions())
    sim.simmd.step(100)


    sim.simmd.context.setParameter('z0',z0)

    return [x0, y0, z0], dz

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

def bondGrapheneSheets(sim, numSheets):
    graph_residues = list(filter(lambda r: r.name == "grph", sim.simmd.topology.residues()))
    print("Restrained atoms ({:d} residues):".format(len(list(graph_residues))))
    restrained_atoms = [[] for y in range(2 * numSheets)]
    for i, res in enumerate(graph_residues):
        for atom in res._atoms:
            if atom.name in map(lambda i: 'C' + str(i), range(0,800,80)):
                restrained_atoms[i].append(atom)

    i = 0
    print(len(restrained_atoms[1]))
    while i < len(restrained_atoms):
        for j in range(len(restrained_atoms[i])):
            # make a bond between this sheet and the adjacent one
            sim.simmd.topology.addBond(restrained_atoms[i][j], restrained_atoms[i + 1][j])
        # skip to the next pair of adjacent sheets
        i += 2

    if len(restrained_atoms) > 2:
        for i in range(len(restrained_atoms[0])):
            # make a bond between the first sheet and the last one
            sim.simmd.topology.addBond(restrained_atoms[0][i], restrained_atoms[-1][i])

    return restrained_atoms
