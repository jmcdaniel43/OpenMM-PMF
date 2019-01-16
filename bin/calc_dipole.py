import mdtraj
from mdtraj import compute_displacements, compute_distances
# for some reason we have to import mdtraj before openmm, not sure why...
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
from numpy import linalg as la
import argparse


def calc_dipole(namesolv, topology, traj, bond_definitions, forcefield, bulk = False, mol_pdb = None):
    """
    This script calculates the distribution of molecular dipole moments.
    It calculates distributions both including and omitting polarization.

    Parameters
    ----------
    namesolv: string
        The name of the molecule to calculate the dipole for.
    topology: string
        Filename of the topology file for the entire trajectory.
    traj: string
        Filename of the trajectory.
    bond_definitions:
        Filename of the forcefield residue definitions.
    forcefield:
        Filename of the forcefield definitions.
    bulk:
        Whether to calculate the molecule dipole as an average over all like residues in the bulk.
    mol_pdb:
        If the full system topology contains multiple forcefield files (or is complex for some othe reason),
        it may be easier just to specify a single pdb file for the molecule whose dipole is being calculated.
        In this case, give the full system topology in the usual argument, and give a supplemental pdb filename
        using this argument to calculate charges on the single molecule.

        This is necessary because the force field files are only used for determining charges on the molecule,
        but after charges are determined the forcefield files are no longer used.
    """

    temperature=300

    # use OpenMM library to get charges
    # here we load pdb without drude oscillators
    if mol_pdb is not None:
        pdb = PDBFile(mol_pdb)
    else:
        pdb = PDBFile(topology)
    pdb.topology.loadBondDefinitions(bond_definitions)
    pdb.topology.createStandardBonds();
    modeller = Modeller(pdb.topology, pdb.positions)

    forcefield = ForceField(forcefield)

    modeller.addExtraParticles(forcefield)
    system = forcefield.createSystem(modeller.topology, nonbondedCutoff=1.4*nanometer, constraints=None, rigidWater=True)
    nbondedForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == NonbondedForce][0]

    # this is not important, this is just to create openmm object that we can use to access topology
    integ_md = DrudeLangevinIntegrator(temperature, 1/picosecond, 1*kelvin, 1/picosecond, 0.001*picoseconds)
    platform = Platform.getPlatformByName('CPU')
    simmd = Simulation(modeller.topology, system, integ_md, platform)

    #####
    #
    # construct charge array for solvent molecule
    # also construct charge array for static charges,
    # so that we can separately compute static and induced dipole moments
    #
    ####
    charges=[]
    chargesnopol=[]
    for res in simmd.topology.residues():
        if res.name == namesolv:
            for i in range(len(res._atoms)):
                index = res._atoms[i].index
                (q, sig, eps) = nbondedForce.getParticleParameters(index)
                charges.append( q._value )
                if 'D' in res._atoms[i].name:
                # no charge on drude without polarization
                    chargesnopol.append(0.0)
                else:
                # to get non-polarization charges,
                # look for matching drude particle
                    typeA=res._atoms[i].name
                    typeD= "D" + typeA
                    indexD=-1
                    for j in range(len(res._atoms)):
                        if typeD == res._atoms[j].name:
                              indexD = res._atoms[j].index
                              break
                    if indexD >= 0:
                        (qD , sig, eps) = nbondedForce.getParticleParameters(indexD)
                        qstatic = q._value + qD._value
                    else:
                        qstatic = q._value
                    chargesnopol.append( qstatic )
            break

    systemTopology = mdtraj.load(topology).topology

    # get indices solvent molecules
    solvent=[]
    for res in systemTopology.residues:
        if res.name == namesolv:
            solvent.append(list(map(lambda a: a.index, res._atoms)))
            if not bulk:
                break

    center_of_charge_atoms = []
    for res in systemTopology.residues:
        if res.name == namesolv:
            for atom in res._atoms:
                atom_name = {
                    "BF4": "B",
                    "TMA": "N"
                }[namesolv]

                if atom.name == atom_name:
                    center_of_charge_atoms.append(atom.index)
                    break

            if not bulk:
                break


    # convert from e*nm to Debye
    conv=1.0/0.0208194

    framestart=0
    frameend=framestart+1

    muhist=[]
    munopolhist=[]
    for i in range(framestart,frameend):
        frame = mdtraj.load_frame(traj, i, top=topology)

        # calculate dipole of each solvent molecule
        for j in range(len(solvent)):
            # make a list of all atom pairs with the charge-carrying center
            res_center = center_of_charge_atoms[j]
            atom_pairs = [[res_center, index] for index in solvent[j]]
            disp = compute_displacements(frame, atom_pairs, periodic = True)

            mu = np.array([0.0, 0.0, 0.0])
            munopol = np.array([0.0, 0.0, 0.0])
            for k in range(len(solvent[j])):
                munopol += chargesnopol[k] * disp[0][k]
                mu += charges[k] * disp[0][k]

            muhist.append(la.norm(mu) * conv)
            munopolhist.append(la.norm(munopol) * conv)


    hist1 = np.histogram(muhist, 50, density=True)
    hist2 = np.histogram(munopolhist, 50, density=True)

    print("dipole distribution with polarization")
    for i in range(len(hist1[0])):
       print (hist1[1][i], hist1[0][i])

    print("dipole distribution without polarization")
    for i in range(len(hist2[0])):
       print (hist2[1][i], hist2[0][i])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("mol", help="residue name of molecule to calculate dipole")
    parser.add_argument("top", help="pdb topology file")
    parser.add_argument("traj", help="trajectory file")
    parser.add_argument("bond_definitions", help="forcefield residue file")
    parser.add_argument("forcefield", help="forcefield definition file")
    parser.add_argument("--bulk", action="store_true", help="calculate ion dipole in the bulk")
    parser.add_argument("--mol_pdb", help="if the total system has complicated forcefield requirements, it may be easier to find the molecular charge information via a separate pdb file")
    args = parser.parse_args()

    calc_dipole(args.mol, args.top, args.traj, args.bond_definitions, args.forcefield, bulk = args.bulk, mol_pdb = args.mol_pdb)

