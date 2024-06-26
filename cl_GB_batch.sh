#!/bin/bash

#SBATCH --job-name=GB_DM
#SBATCH --time=1-0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --partition=hpc
#SBATCH --output=/home/tomas_mendes/projects/cosmolattice/cosmo.out
#SBATCH --error=/home/tomas_mendes/projects/cosmolattice/cosmo.err

RUNPATH=/home/tomas_mendes/projects/cosmolattice
cd $RUNPATH

#WHAT IS THIS?
module purge

module load gcc-11.3
module load gcc11/openmpi/4.1.4
module load gcc11/libs/fftw/3.3.10

n_lattice=$1
t_max=$2
Hkin_val=$3
cutoff=$4
m_val=$5
k_min=$6
k_max=$7
n_run=$8
dir_index=${9}

srun --mpi=openmpi ./GB_wMass outputfile=./run_${dir_index}/ expansion=true evolver=VV4 fixedBackground=true omegaEoS=0.333 H0=${Hkin_val} N=${n_lattice} dt=0.1 kIR=${k_min} tOutputFreq=0.1 tOutputInfreq=1000 tOutputRareFreq=1000 tOutputVerb=10 tMax=${t_max} PS_type=2 PS_version=1 GWprojectorType=1 withGWs=false kCutOff=${k_max} initial_amplitudes=0 initial_momenta=0 M=${m_val} Hkin=${Hkin_val} cutoff=${cutoff} baseSeed=1922




