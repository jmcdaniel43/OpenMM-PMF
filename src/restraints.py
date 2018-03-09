from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

def findCenterOfPore(sim):
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
    state = sim.simmd.context.getState(getPositions=True)
    prePotentialPositions = state.getPositions()

    for i, center in enumerate(atomIndices):
        for index in center:
            coords[i][0] += prePotentialPositions[index][0] / nanometer
            coords[i][1] += prePotentialPositions[index][1] / nanometer
            coords[i][2] += prePotentialPositions[index][2] / nanometer

        coords[i][0] /= len(center)
        coords[i][1] /= len(center)
        coords[i][2] /= len(center)

    return coords[0], coords[1], coords[2]

def addIonUmbrellaPotential(sim, distanceFromPore, numbrella, kxy):
    """
    Apply an umbrella potential window at the specified distanceFromPore
    from the edge of the pore for the Boron atom in the first BF4 residue.

    Parameters
    _________
    distanceFromPore: float
        Distance from the edge of the pore for the potential to exist.
    numbrella: int
        Number of umbrella potentials to be used
    kxy: float
        Force constant for the potential in the x, y dimensions.

    Returns
    ______

    potentialCoordinates, [x0,y0,z0]: [float]
        Starting location of the potential
    dz: float
        Distance to move the umbrella after each iteration.
    """

    print('Setting umbrella forces about equilibrated pore')

    # Locate ion
    name='BF4'
    aname='B'
    index=-1
    for res in sim.simmd.topology.residues():
        if res.name == name:
            index = res._atoms[0].index
            break
    print("Ion index:", index)

    centerCoords, interiorCoords, exteriorCoords = findCenterOfPore(sim)
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
    ZForce.addGlobalParameter('z0', z0)
    ZForce.addGlobalParameter('kz', kz)

    # PMF apply umbrella potential to x, y directions
    XYForce = CustomExternalForce("0.5*kxy*periodicdistance(x,y,z,x0,y0,z)^2")
    sim.system.addForce(XYForce)
    XYForce.addParticle(index)
    XYForce.addGlobalParameter('x0', x0)
    XYForce.addGlobalParameter('y0', y0)
    XYForce.addGlobalParameter('kxy', kxy)

    return index, [x0, y0, z0], dz

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
