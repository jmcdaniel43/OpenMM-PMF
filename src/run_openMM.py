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

def main(filename = "md_nvt_prod"):
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="brief description of the simulation")
    parser.add_argument("pdb", help="PDB file with initial positions")
    parser.add_argument("sheets", type=int, help="number of graphene sheets per electrode", choices=[1,2])
    parser.add_argument("pore", type=int, help="pore size", choices=[7,10,14])
    parser.add_argument("--ffdir", default="ffdir", help="directory containing force field files")
    parser.add_argument("--kxy", type=float, default=5, help="x-y dimension spring constant for ion constraint")
    parser.add_argument("--prefix", default="md_nvt_prod", help="prefix for output files (*.dcd, *.pdb)")
    parser.add_argument("--eqNPT", default=5000, type=int, help="how long to equilibrate NPT")
    parser.add_argument("--eqNVT", default=5000, type=int, help="how long to equilibrate NVT")
    parser.add_argument("--longSample", default="false", help="whether to sample 120ns of 0.02 nm windows")
    args = parser.parse_args()
    print(args)


    outdir = "{:d}sheet/{:d}pore/output_{:s}/".format(args.sheets, args.pore, args.name)
    os.makedirs(outdir)
    print(outdir)

    ffdir = args.ffdir + "/"
    sim = PMF_Simulation(
        args.pdb,
        filename,
        outdir,
        bondDefinitionFiles = [
            ffdir+'sapt_residues.xml',
            ffdir+'graph_residue_{:d}.xml'.format(args.pore)
        ], forceFieldFiles = [
            ffdir+'sapt.xml',
            ffdir+'graph_{:d}.xml'.format(args.pore)
        ]
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

    distanceFromPore = 1.0
    numbrella=60

    if args.longSample.lower() in ["t", "tr", "tru", "true"]:
        print("Running 120ns of 0.02 nm windows")
        numbrella=120

    index, poreCenter, dz = addIonUmbrellaPotential(sim, distanceFromPore, numbrella, args.kxy)
    print("dz", dz)
    z0 = poreCenter[2]

    restrainGrapheneSheets(sim, restrained_atoms, z0 + distanceFromPore)

    boxDims = sim.simmd.topology.getUnitCellDimensions()
    print("Box dimensions: ", boxDims)
    if z0 < 0:
        print("First umbrella is outside periodic box. Moving z coordinate to other side of image.")
        z0 = boxDims[2] / nanometer + z0
    print("Center of umbrella potential (nm)", poreCenter)

    state = sim.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True)
    sim.simmd.context.reinitialize()
    sim.simmd.context.setPositions(state.getPositions())

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
       for i in range(10000):
           # print position of Boron Atom; first print is starting position
           state = sim.simmd.context.getState(getPositions=True, enforcePeriodicBox=True)
           position = state.getPositions()
           print(position[index])

           sim.simmd.step(100)

       z0 = z0 + dz
       sim.simmd.context.setParameter('z0',z0)


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
