from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

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

class PMF_Simulation:
    def __init__(self, pdb_name, filename, outdir, bondDefinitionFiles, forceFieldFiles, gpuDeviceIndex = '0', trajFreq=50000):
        self.temperature = 300 * kelvin
        self.outdir = outdir
        self.filename = filename
        pdb = PDBFile(pdb_name)

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


        platform = Platform.getPlatformByName('CUDA')
        print("Using GPU device:", gpuDeviceIndex)
        #properties = {'OpenCLPrecision': 'mixed', 'OpenCLDeviceIndex': gpuDeviceIndex}
        properties = {'Precision': 'mixed'}

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
            DCDReporter(self.outdir + self.filename + '.dcd', trajFreq),
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

        state = self.simmd.context.getState(getPositions=True, enforcePeriodicBox=True)
        boxVecs = state.getPeriodicBoxVectors()
        self.system.setDefaultPeriodicBoxVectors(boxVecs[0], boxVecs[1], boxVecs[2])
        self.simmd.topology.setPeriodicBoxVectors(boxVecs)
        PDBFile.writeFile(self.simmd.topology, state.getPositions(), open(self.outdir + self.filename + '_equil_npt.pdb', 'w'))

    def equilibrate_nvt(self, time_ps = 5000, checkpointFile = None):
        print('Starting NVT Equilibration ({:d} ps)...'.format(time_ps))

        if checkpointFile is not None:
            print('Using checkpoint, ' + checkpointFile)
            self.simmd.loadCheckpoint(checkpointFile)
            return

        for i in range(time_ps):
            self.simmd.step(1000)

        state = self.simmd.context.getState(getPositions=True, enforcePeriodicBox=True)
        PDBFile.writeFile(self.simmd.topology, state.getPositions(), open(self.outdir + self.filename + '_equil_nvt.pdb', 'w'))

