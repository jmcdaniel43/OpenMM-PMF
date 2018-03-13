# OpenMM-PMF
Code for computing PMFs.  Directly used to compute ion diffusion barriers through pores in graphene sheets.

Creating systems
-----------

Electrolyte solution pdb's should be created in the `setup_data` directory. Create the pdb with packmol using one of the input scripts available. In order to work with the next script, these systems should be labeled `<system_name>_<solvent_name>`.

Then in the top-level directory, use `system_setup.sh` with two arguments (solvent and system) to create a collection of supercapacitor pdb files in the `pdb` directory. This script combines the electrolyte solution with eelectrode pdb's to create the production system.

```bash
$ cd system_data
$ packmol < bf4_tma_tmea_acn.inp # bf4_tma_tmea_acn.pdb is created
cd ..

./system_setup.sh acn bf4_tma_tmea # pdb/acn/bf4_tma_tmea is created with a collection of pore sizes and electrode sheet numbers
```
Use with SLURM queueing systems
--------------

`make_slurm.sh` makes the set of slurm scripts for all combinations of solvent, system, ion to diffuse, and electrode pore size/sheet count.

Running the simulation
----------

The simulation can be run with
```bash
$ src/run_openMM.py bf4_tma_tmea acn bf4 1 10 |& tee output_logs/bf4_tma_tmea_acn_1_10_bf4diff.log
```
This will create necessary directories to produce output in a directory `<sheet count>/<solvent>/<pore size>/output_<system>_<diffusing ion>diff[_<tag>]`.

Since there is a significant amount of output, it is probably a good idea to pipe it to another file for replication and analysis.

---

Electrolyte solutions
-----------

The available solution components are listed if you run the `--help` switch on `src/run_openMM.py`. These are the current availabe components.

__Electrolytes__

- bf4
- tma
- tmea

__Solvents__

- acn
- dce

---

Analysis
-------

## PMF

For analysis, there is a script to calculate PMF of the ion diffusion. It can be run like this:
```bash
$ mkdir pmf
$ cd pmf
$ convert_umbrella_output.py
USAGE: <input file> <spring constant> <dz> <nstep> <numbrella>
$ convert_umbrella_output.py 2_10_dce.log 2000 0.4 10000 60 > whaminput
# then the pmf can be calculated with
$ wham 62.0 86.0 15 0.01 300.0 0 whaminput pmf 1000 143289

```
- `input file`: the output of the production script
- `spring constant`: this should be 2000 unless the production script was changed.
- `dz`: the distance the umbrella potential moved between each window 
- `nstep`: the number of location printouts per window from the script
- `numbrella`: number of potential windows utilized over the course of the simulation

## Solvation

The solvation of the ion through the pore is measured (per window) by the coordination number found by integration of the radial distribution function.

The script to do this is as follows:
```bash
rdf.py traj_dce --bulk=t | rdf_smoother.py
```

The first argument is the directory containing the tracjetory from the simulation. A switch is provided to calculate the bulk coordination number (omit to only calculate for the diffusion ion). The output of this is piped to a smoothing function that cleans up noise from the original rdf.

