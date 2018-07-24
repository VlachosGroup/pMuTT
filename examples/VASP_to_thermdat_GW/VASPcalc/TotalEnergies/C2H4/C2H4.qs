#!/bin/bash
#$ -cwd
#$ -pe mpi 20
#$ -l exclusive=1
#$ -o C2H4.stdout
#$ -e C2H4.stderr
#$ -l m_mem_free=2G
#$ -l h_cpu=36:00:00
#$ -m beas
#$ -M ggu@udel.edu

source /etc/profile.d/valet.sh
vpkg_rollback all
vpkg_require vasp/5.3.2:d3,gamma_only,intel,mpi,vtst


mpiexec -n 20 /home/work/ccei_biomass/sw/vasp/5.3.2-intel64-openmpi+VTST+D3+GAMMA/vasp >& C2H4.out
