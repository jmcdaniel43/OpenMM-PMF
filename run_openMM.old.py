#!/usr/bin/env python

from __future__ import print_function
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from time import gmtime, strftime
from datetime import datetime
import sys
import argparse

class PMF_Simulation:
    def __init__(self, pdb, filename, outdir, bondDefinitionFiles, forceFieldFiles):
        self.temperature = 300 * kelvin
        self.outdir = outdir
        self.filename = filename

        integ_md = DrudeLangevinIntegrator(self.temperature, 1/picosecond, 1*kelvin, 1/picosecond, 0.001*picoseconds)
        integ_md.setMaxDrudeDistance(0.02)  # this should prevent polarization catastrophe during equilibration, but shouldn't affect results afterwards ( 0.2 Angstrom displacement is very large for equil. Drudes)

        for ff in bondDefinitionFiles:
            pdb.topology.loadBondDefinitions(ff)
        pdb.topology.createStandardBonds();

        modeller = Modeller(pdb.topology, pdb.positions)
        forcefield = ForceField(*forceFieldFiles)
        modeller.addExtraParticles(forcefield)

        self.system = forcefield.createSystem(modeller.topology, nonbondedCutoff=1.4*nanometer, constraints=None, rigidWater=True)
        self.nbondedForce = [f for f in [self.system.getForce(i) for i in range(self.system.getNumForces())] if type(f) == NonbondedForce][0]
        self.customNonbondedForce = [f for f in [self.system.getForce(i) for i in range(self.system.getNumForces())] if type(f) == CustomNonbondedForce][0]
        self.custombond = [f for f in [self.system.getForce(i) for i in range(self.system.getNumForces())] if type(f) == CustomBondForce][0]
        self.drudeForce = [f for f in [self.system.getForce(i) for i in range(self.system.getNumForces())] if type(f) == DrudeForce][0]
        self.nbondedForce.setNonbondedMethod(NonbondedForce.PME)
        self.customNonbondedForce.setNonbondedMethod(min(self.nbondedForce.getNonbondedMethod(),NonbondedForce.CutoffPeriodic))
        self.customNonbondedForce.setUseLongRangeCorrection(True)

        print('nbMethod : ', self.customNonbondedForce.getNonbondedMethod())

        for i in range(self.system.getNumForces()):
            f = self.system.getForce(i)
            type(f)
            f.setForceGroup(i)
            # Here we are adding periodic boundaries to intra-molecular interactions.  Note that DrudeForce does not have this attribute, and
            # so if we want to use thole screening for graphite sheets we might have to implement periodic boundaries for this force type
            if type(f) == HarmonicBondForce or type(f) == HarmonicAngleForce or type(f) == PeriodicTorsionForce or type(f) == RBTorsionForce:
               f.setUsesPeriodicBoundaryConditions(True)
            f.usesPeriodicBoundaryConditions()

        totmass = 0.*dalton
        for i in range(self.system.getNumParticles()):
            totmass += self.system.getParticleMass(i)


        platform = Platform.getPlatformByName('OpenCL')
        properties = {'OpenCLPrecision': 'mixed','OpenCLDeviceIndex':'0'}

        self.simmd = Simulation(modeller.topology, self.system, integ_md, platform, properties)
        self.simmd.context.setPositions(modeller.positions)

        platform = self.simmd.context.getPlatform()
        platformname = platform.getName();
        print(platformname)

        # Print initial energy
        state = self.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True)
        print('Initial Energy')
        print(str(state.getKineticEnergy()))
        print(str(state.getPotentialEnergy()))

        for i in range(self.system.getNumForces()):
            f = self.system.getForce(i)
            print(type(f), str(self.simmd.context.getState(getEnergy=True, groups=2**i).getPotentialEnergy()))
            f.setForceGroup(i)

        # write initial pdb with drude oscillators
        position = state.getPositions()
        self.simmd.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
        PDBFile.writeFile(self.simmd.topology, position, open(self.outdir + self.filename + '_start_drudes.pdb', 'w'))

        self.simmd.reporters = [
            DCDReporter(self.outdir + self.filename + '.dcd', 50000),
            CheckpointReporter(self.outdir + self.filename + '.chk', 10000)
        ]
        self.simmd.reporters[1].report(self.simmd,state)

    def equilibrate_npt(self, time_ps = 5000, checkpointFile = None):
        print('Starting NPT Equilibration ({:d} ps)...'.format(time_ps))

        if checkpointFile is not None:
            print('Using checkpoint, ' + checkpointFile)
            self.simmd.loadCheckpoint(checkpointFile)
            return

        pressure = Vec3(1.0,1.0,1.0)*atmosphere
        barofreq = 100

        # allow only the z-dimension to change, graphene x/y layer is fixed
        barostat = MonteCarloAnisotropicBarostat(pressure,self.temperature,False,False,True,barofreq)
        barostatForceIndex = self.system.addForce(barostat)

        state = self.simmd.context.getState(getPositions=True)
        initialPositions = state.getPositions()
        self.simmd.context.reinitialize()
        self.simmd.context.setPositions(initialPositions)

        barofreq = barostat.getFrequency()
        print("barofreq: ", barofreq)

        for i in range(time_ps):
            self.simmd.step(1000)


        self.system.removeForce(barostatForceIndex)

        state = self.simmd.context.getState(getPositions=True)
        postNPTPositions = state.getPositions()
        self.simmd.context.reinitialize()
        self.simmd.context.setPositions(postNPTPositions)
        PDBFile.writeFile(self.simmd.topology, postNPTPositions, open(self.outdir + self.filename + '_equil_npt.pdb', 'w'))

    def equilibrate_nvt(self, time_ps = 5000, checkpointFile = None):
        print('Starting NVT Equilibration ({:d} ps)...'.format(time_ps))

        if checkpointFile is not None:
            print('Using checkpoint, ' + checkpointFile)
            self.simmd.loadCheckpoint(checkpointFile)
            return

        for i in range(time_ps):
            self.simmd.step(1000)

        state = self.simmd.context.getState(getPositions=True)
        position = state.getPositions()
        self.simmd.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
        PDBFile.writeFile(self.simmd.topology, position, open(self.outdir + self.filename + '_equil_nvt.pdb', 'w'))

    def addIonUmbrellaPotential(self, kxy):
        print('Setting umbrella forces about equilibrated pore')

        ####### ################
        # Find center of pore
        ####### ###############

        # Find relevant C atoms on the edge of the pore
        poreEdgeAtomIndicies = []
        sheetResIndex = -1
        print("Pore edge atom indicies:", end=" ")
        for res in self.simmd.topology.residues():
            if res.name == "grph":
                if sheetResIndex == -1:
                    # we only want to consider atoms from one sheet of the graphene
                    # so force one residue index for this loop
                    sheetResIndex = res.index
                elif res.index != sheetResIndex:
                    continue

                for atom in res._atoms:
                    if atom.name in ["C354", "C274", "C148", "C68"]:
                        print(atom.index, end=" ")
                        poreEdgeAtomIndicies.append(atom.index)
        print()

        state = self.simmd.context.getState(getPositions=True)
        prePotentialPositions = state.getPositions()
        x0 = 0
        y0 = 0
        z0 = 0

        for index in poreEdgeAtomIndicies:
            x0 += prePotentialPositions[index][0] / nanometer
            y0 += prePotentialPositions[index][1] / nanometer
            z0 += prePotentialPositions[index][2] / nanometer

        x0 /= len(poreEdgeAtomIndicies)
        y0 /= len(poreEdgeAtomIndicies)
        z0 /= len(poreEdgeAtomIndicies)

        ####### ################
        # Locate ion
        ####### ################

        name='BF4'
        aname='B'
        index=-1
        for res in self.simmd.topology.residues():
            if res.name == name:
                index = res._atoms[0].index
                break

        print("Ion index:", index)

        ####### ################
        # Apply potential
        ####### ################

        # PMF apply umbrella potential to z direction
        kz=2000.0  # kJ/mol/nm^2
        ZForce = CustomExternalForce("0.5*kz*periodicdistance(x,y,z,x,y,z0)^2")
        self.system.addForce(ZForce)
        ZForce.addParticle(index)
        ZForce.addGlobalParameter('z0', z0)
        ZForce.addGlobalParameter('kz', kz)

        # PMF apply umbrella potential to x, y directions
        XYForce = CustomExternalForce("0.5*kxy*periodicdistance(x,y,z,x0,y0,z)^2")
        self.system.addForce(XYForce)
        XYForce.addParticle(index)
        XYForce.addGlobalParameter('x0', x0)
        XYForce.addGlobalParameter('y0', y0)
        XYForce.addGlobalParameter('kxy', kxy)

        return index, x0, y0, z0

    def restrainGrapheneSheet(self, restrained_atoms, z0graph):
        state = self.simmd.context.getState(getPositions=True)
        preGrapheneRestraintPositions = state.getPositions()

        all_atoms = [item for sublist in restrained_atoms for item in sublist]
        all_atoms = map(lambda x: x.index, all_atoms)
        for i, carbon in enumerate(all_atoms):
            x0graph = preGrapheneRestraintPositions[carbon][0] / nanometer
            y0graph = preGrapheneRestraintPositions[carbon][1] / nanometer

            # add restraining force to 10 atoms in the sheet
            kgraph=2000.0  # kJ/mol/nm^2
            GraphForce = CustomExternalForce("0.5*kgraph*periodicdistance(x,y,z, x0graph" + str(i) + ", y0graph" + str(i) + ", z0graph)^2")
            self.system.addForce(GraphForce)
            GraphForce.addParticle(carbon)
            GraphForce.addGlobalParameter('kgraph', kgraph)
            GraphForce.addGlobalParameter('x0graph' + str(i), x0graph)
            GraphForce.addGlobalParameter('y0graph' + str(i), y0graph)
            GraphForce.addGlobalParameter('z0graph', z0graph)

def main(filename = "md_nvt_prod"):
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="brief description of the simulation")
    parser.add_argument("pdb", help="PDB file with initial positions")
    parser.add_argument("pore", type=int, help="pore size", choices=[7,10,14])
    parser.add_argument("sheets", type=int, help="number of graphene sheets per electrode", choices=[1,2])
    parser.add_argument("--ffdir", default="ffdir", help="directory containing force field files")
    parser.add_argument("--kxy", type=float, default=5, help="x-y dimension spring constant for ion constraint")
    parser.add_argument("--prefix", default="md_nvt_prod", help="prefix for output files (*.dcd, *.pdb)")
    parser.add_argument("--eqNPT", default="true", help="whether to equilibrate NPT")
    parser.add_argument("--eqNVT", default="true", help="whether to equilibrate NVT")
    args = parser.parse_args()
    print(args)


    outdir = "{:d}sheet/{:d}pore/output_{:s}/".format(args.sheets, args.pore, args.name)
    os.makedirs(outdir)
    print(outdir)

    ffdir = args.ffdir + "/"
    sim = PMF_Simulation(
        PDBFile(args.pdb),
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

    if args.eqNPT.lower() in ["t", "tr", "tru", "true"]:
        #sim.equilibrate_npt(checkpointFile = "md_npt_equil_bonded.chk")
        sim.equilibrate_npt(5000)
    else:
        print("Skipping NPT equilibration")

    state = sim.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True)
    postNPTPositions = state.getPositions()
    print('Post-NPT-equilibration Energy')
    print(str(state.getKineticEnergy()))
    print(str(state.getPotentialEnergy()))

    for i in range(sim.system.getNumForces()):
        f = sim.system.getForce(i)
        print(type(f), str(sim.simmd.context.getState(getEnergy=True, groups=2**i).getPotentialEnergy()))
        f.setForceGroup(i)


    graph_residues = list(filter(lambda r: r.name == "grph", sim.simmd.topology.residues()))
    print("Restrained atoms ({:d} residues):".format(len(list(graph_residues))))
    restrained_atoms = [[]] * 2 * args.sheets
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

    sim.simmd.context.reinitialize()
    sim.simmd.context.setPositions(postNPTPositions)

    if args.eqNVT.lower() in ["t", "tr", "tru", "true"]:
        #sim.equilibrate_nvt(500, checkpointFile = "equil/equil_npt.chk")
        sim.equilibrate_nvt(5000)
    else:
        print("Skipping NVT equilibration")

    index, x0, y0, z0 = sim.addIonUmbrellaPotential(args.kxy)

    z0graph = z0
    z0 = z0 - 1 # to start 10 A away from the pore

    boxDims = sim.simmd.topology.getUnitCellDimensions()
    print("Box dimensions: ", boxDims)
    if z0 < 1:
        z0 = boxDims[2] / nanometer + (z0 - 1)
    print("Center of umbrella potential (nm): x0 y0 z0", x0, y0, z0)

    sim.restrainGrapheneSheet(restrained_atoms, z0graph)

    state = sim.simmd.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True)
    postNVTPositions = state.getPositions()
    print('Post-NVT-equilibration Energy')
    print(str(state.getKineticEnergy()))
    print(str(state.getPotentialEnergy()))

    for i in range(sim.system.getNumForces()):
        if i >= 31:
            print("Forces after this cannot be printed with this function: {:d} forces not printed".format(sim.system.getNumForces() - i))
            break
        f = sim.system.getForce(i)
        print(type(f), str(sim.simmd.context.getState(getEnergy=True, groups={i}).getPotentialEnergy()))
        f.setForceGroup(i)

    sim.simmd.context.reinitialize()
    sim.simmd.context.setPositions(postNVTPositions)

    t1 = datetime.now()
    print('Starting Production NVT PMF Simulation...')
    # loop over umbrella positions
    dz=0.04
    numbrella=60
    for iu in range(numbrella):
       z0 = z0 + dz
       sim.simmd.context.setParameter('z0',z0)

       for i in range(10000):
           sim.simmd.step(100)

           # print position of Boron Atom
           state = sim.simmd.context.getState(getPositions=True)
           position = state.getPositions()
           print(position[index])


    t2 = datetime.now()
    t3 = t2 - t1
    print('simulation took', t3.seconds,'seconds')
    state = sim.simmd.context.getState(getEnergy=True,getForces=True,getPositions=True,enforcePeriodicBox=True)
    position = state.getPositions()
    sim.simmd.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
    PDBFile.writeFile(sim.simmd.topology, position, open(outdir + filename + '.pdb', 'w'))

    print('Done!')

if __name__ == '__main__':
    main()
