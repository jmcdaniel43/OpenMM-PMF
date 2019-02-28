#!/usr/bin/env python

from __future__ import print_function
from sys import stdout
from time import gmtime, strftime
from datetime import datetime
import sys
import argparse
import os
import re

# from simtk.openmm.app.vvintegrator import Barostat
from simtk.unit import *

from simulation import PMF_Simulation, cations, anions, solvents
from restraints import *

default_eqNPT = 5000
default_eqNVT = 5000

def main(args, filename = "md_nvt_prod"):
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

    pdb = "pdb/{:s}/{:s}/SC_start_{:d}_{:d}.pdb".format(args.solvent, args.system, args.sheets, args.pore)
    if args.pdb is not None:
        if not os.path.exists(args.pdb):
            print("pdb files does not exist ({:s})".format(args.pdb))
            print(dir(args.eqNVT))
            print("assuming equilibration should be reset to default: eqNPT={:d} eqNVT={:d}".format(default_eqNPT, default_eqNVT))
            args.eqNPT = default_eqNPT
            args.eqNVT = default_eqNVT
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

    state = sim.simmd.context.getState(getPositions=True,enforcePeriodicBox=True)
    print('Post-NPT equilibration Energy')
    printEnergy(sim)

    restrained_atoms = bondGrapheneSheets(sim, args.sheets)

    sim.simmd.context.reinitialize()
    sim.simmd.context.setPositions(state.getPositions())

    if args.eqNVT > 0:
        sim.equilibrate_nvt(args.eqNVT)
    else:
        print("Skipping NVT equilibration")

    if args.eqOnly:
        return

    distanceFromPore = args.startingDistance
    numbrella=60

    if args.longSample:
        print("Running 120ns of 0.02 nm windows")
        numbrella=120

    index = None
    if args.restart_ion_index:
        index = args.restart_ion_index

    potentialCenter, dz, index = addIonUmbrellaPotential(sim, distanceFromPore, numbrella, args.kxy, args.ion, index = index)
    print("dz", dz)
    z0 = potentialCenter[2]


    restrainGrapheneSheets(sim, restrained_atoms, z0 + distanceFromPore)

    boxDims = sim.simmd.topology.getUnitCellDimensions()
    print("Box dimensions: ", boxDims)
    print("Center of umbrella potential (x,y,z (nm)):", potentialCenter[0], potentialCenter[1], potentialCenter[2])

    state = sim.simmd.context.getState(getEnergy=True, getForces=True, getVelocities=True, getPositions=True, enforcePeriodicBox=True)
    sim.simmd.context.reinitialize()
    sim.simmd.context.setPositions(state.getPositions())

    if args.restart_window > 0:
        z0 += args.restart_window*dz
        sim.simmd.context.setParameter('z0',z0)
        print("restart window: moving window z0 to", z0, ", window", args.restart_window)

    # equilibrate 100 picoseconds to ensure ion is in the correct umbrella
    print("Equilibrating to first umbrella window")
    sim.simmd.step(100000)

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
    for iu in range(args.restart_window, numbrella):
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
    state = sim.simmd.context.getState(getEnergy=True, getForces=True, getPositions=True, enforcePeriodicBox=True)
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
    parser.add_argument("--eqNPT", default=default_eqNPT, type=int, help="how long to equilibrate NPT")
    parser.add_argument("--eqNVT", default=default_eqNVT, type=int, help="how long to equilibrate NVT")
    parser.add_argument("--longSample", action="store_true", help="whether to sample 120ns of 0.02 nm windows")
    parser.add_argument("--gpuDevice", default='0', type=str, help="gpu device index")
    parser.add_argument("--trajFreq", default=50000, type=int, help="how often trajectory frames are printed (fs)")
    parser.add_argument("-f", "--forceOverwrite", action="store_true", help="force overwrite of existing data")
    parser.add_argument("--windowSize", default=10000, type=int, help="length of time to spend in each potential window (units of 0.1 ps)")
    parser.add_argument("--eqOnly", action="store_true", help="use if you want to end the simulation after equilibration")
    parser.add_argument("--startingDistance", type=float, default=1.0, help="distance from the edge of the pore to start the umbrella sampling")
    parser.add_argument("--freezeSheets", action="store_true", help="freeze graphene sheets in place after equilibration")
    parser.add_argument("--eqExcludeGraphene", action="store_true", help="exclude graphene sheets in NPT equilibration")
    parser.add_argument("--restart_window", type=int, default=0, help="restart the trajectory from a specific umbrella window (zero-indexed)")
    parser.add_argument("--restart_ion_index", type=int, default=1, help="restart the trajectory with a specific ion index (you probably want to do this)")

    args = parser.parse_args()
    print(args)

    main(args)
