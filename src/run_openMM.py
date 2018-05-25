#!/usr/bin/env python

from __future__ import print_function
from sys import stdout
from time import gmtime, strftime
from datetime import datetime
import sys
import argparse
import os

from simulation import PMF_Simulation
from restraints import *

solvents = ["acn", "dce"]
anions = {
    "BF4": "B",
}
cations = {
    "TMA": "N",
    "TMEA": "N"
}
ionAtomNameMap = {}
ionAtomNameMap.update(anions)
ionAtomNameMap.update(cations)

def main(filename = "md_nvt_prod"):
    parser = argparse.ArgumentParser()
    parser.add_argument("system", help="system of anions and cations to use. Syntax is <anion>_,cation1>[_<cation2>...]")
    parser.add_argument("solvent", help="solvent to use for system", choices=solvents)
    parser.add_argument("ion", help="residue name of ion to apply umbrella potential to", choices=list(ionAtomNameMap.keys()))
    parser.add_argument("sheets", type=int, help="number of graphene sheets per electrode", choices=[1,2])
    parser.add_argument("pore", type=int, help="pore size", choices=[7,10,14])
    parser.add_argument("--tag", help="brief description of the simulation")
    parser.add_argument("--pdb", help="PDB file with initial positions. Default is pdb/<solvent>/<system>/SC_start_<sheets>_<pore>.pdb")
    parser.add_argument("--ffdir", default="ffdir", help="directory containing force field files")
    parser.add_argument("--kxy", type=float, default=5, help="x-y dimension spring constant for ion constraint")
    parser.add_argument("--prefix", default="md_nvt_prod", help="prefix for output files (*.dcd, *.pdb)")
    parser.add_argument("--eqNPT", default=5000, type=int, help="how long to equilibrate NPT")
    parser.add_argument("--eqNVT", default=5000, type=int, help="how long to equilibrate NVT")
    parser.add_argument("--longSample", action="store_true", help="whether to sample 120ns of 0.02 nm windows")
    parser.add_argument("--gpuDevice", default='0', type=str, help="gpu device index")
    parser.add_argument("--trajFreq", default=50000, type=int, help="how often trajectory frames are printed (fs)")
    parser.add_argument("-f", "--forceOverwrite", action="store_true", help="force overwrite of existing data")
    parser.add_argument("--windowSize", default=10000, type=int, help="length of time to spend in each potential window (units of 0.1 ps)")
    parser.add_argument("--eqOnly", action="store_true", help="use if you want to end the simulation after equilibration")

    args = parser.parse_args()
    print(args)


    outdir = "{:d}sheet/{:s}/{:d}pore/output_{:s}_{:s}diff".format(args.sheets, args.solvent, args.pore, args.system, args.ion)
    if args.tag is not None:
        outdir += "_" + args.tag
    outdir += "/"
    try:
        os.makedirs(outdir)
    except FileExistsError as e:
        if not args.forceOverwrite:
            raise e
    print(outdir)

    if args.pdb is None:
        pdb = "pdb/{:s}/{:s}/SC_start_{:d}_{:d}.pdb".format(args.solvent, args.system, args.sheets, args.pore)
    else:
        pdb = args.pdb
    print("Using pdb: ", pdb)

    ffdir = args.ffdir + "/"
    sim = PMF_Simulation(
        pdb,
        filename,
        outdir,
        bondDefinitionFiles = [
            ffdir+'sapt_residues.xml',
            ffdir+'graph_residue_{:d}.xml'.format(args.pore)
        ],
        forceFieldFiles = [
            ffdir+'sapt.xml',
            ffdir+'graph_{:d}.xml'.format(args.pore)
        ],
        gpuDeviceIndex = args.gpuDevice,
        trajFreq = args.trajFreq
    )

    if args.eqNPT > 0:
        sim.equilibrate_npt(args.eqNPT)
    else:
        print("Skipping NPT equilibration")

    state = sim.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True,enforcePeriodicBox=True)
    print('Post-NPT equilibration Energy')
    printEnergy(sim)

    for i in range(sim.system.getNumForces()):
        f = sim.system.getForce(i)
        f.setForceGroup(i)
        print(type(f), str(sim.simmd.context.getState(getEnergy=True, groups={i}).getPotentialEnergy()))

    restrained_atoms = bondGrapheneSheets(sim, args.sheets)

    sim.simmd.context.reinitialize()
    sim.simmd.context.setPositions(state.getPositions())

    if args.eqNVT > 0:
        sim.equilibrate_nvt(args.eqNVT)
    else:
        print("Skipping NVT equilibration")

    if args.eqOnly:
        return

    distanceFromPore = 1.0
    numbrella=60

    if args.longSample:
        print("Running 120ns of 0.02 nm windows")
        numbrella=120

    # Locate ion
    aname = ionAtomNameMap[args.ion]
    index=-1
    for res in sim.simmd.topology.residues():
        if res.name == args.ion:
            for atom in res._atoms:
                if atom.name == aname:
                    print(str(atom))
                    index = atom.index
                    break
            break
    print("Ion index:", index)
    if index == -1:
        raise IndexError("ion index was not found in the system")

    poreCenter, dz = addIonUmbrellaPotential(sim, distanceFromPore, numbrella, args.kxy, index, args.pore)
    print("dz", dz)
    z0 = poreCenter[2]

    restrainGrapheneSheets(sim, restrained_atoms, z0 + distanceFromPore)

    boxDims = sim.simmd.topology.getUnitCellDimensions()
    print("Box dimensions: ", boxDims)
    print("Center of umbrella potential (nm)", poreCenter)

    state = sim.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True, enforcePeriodicBox=True)
    sim.simmd.context.reinitialize()
    sim.simmd.context.setPositions(state.getPositions())

    # after reinitialization, the parameter gets reset to default value
    # set this again to make make sure it sticks
    sim.simmd.context.setParameter('z0',z0)

    print('Post-NVT (after reinitialize)')
    print('Box vectors:', sim.simmd.context.getState().getPeriodicBoxVectors())
    print("Equilibration Energy")
    printEnergy(sim)


    ###############

    # Start sampling

    ################

    print('Starting Production NVT PMF Simulation...')
    t1 = datetime.now()

    # loop over umbrella positions
    for iu in range(numbrella):
        for i in range(args.windowSize):
            # print position of Boron Atom; first print is starting position
            state = sim.simmd.context.getState(getPositions=True, enforcePeriodicBox=True)
            position = state.getPositions()
            print(position[index])

            sim.simmd.step(100)

        z0 = z0 + dz
        sim.simmd.context.setParameter('z0',z0)
        print("moving window z0 to", z0)


    t2 = datetime.now()
    t3 = t2 - t1
    print('simulation took', t3.seconds,'seconds')
    state = sim.simmd.context.getState(getEnergy=True,getForces=True,getPositions=True,enforcePeriodicBox=True)
    position = state.getPositions()
    sim.simmd.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
    PDBFile.writeFile(sim.simmd.topology, position, open(outdir + filename + '.pdb', 'w'))

    print('Done!')

def printEnergy(sim):
    state = sim.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True)
    print(str(state.getKineticEnergy()))
    print(str(state.getPotentialEnergy()))

    seenForceGroups = []
    for i in range(sim.system.getNumForces()):
        if i >= 31:
            print("Forces after this cannot be printed with this function: {:d} forces not printed".format(sim.system.getNumForces() - i))
            break
        f = sim.system.getForce(i)
        if f.getForceGroup() not in seenForceGroups:
            seenForceGroups.append(f.getForceGroup)
            print(type(f), str(sim.simmd.context.getState(getEnergy=True, groups={f.getForceGroup()}).getPotentialEnergy()))


if __name__ == '__main__':
    main()
