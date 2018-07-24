#!/bin/sh
#
#PBS -N /user/home/gent/vsc400/vsc40063/data/glen/gas-phase/H2
#PBS -o vasp.out
#PBS -e vasp.err
#PBS -q default
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=8

newgrp - gvasp
module load VASP/5.3.3-ictce-4.1.13-vtst-3.0c-20130327-mt-gamma

ORIGDIR=/user/home/gent/vsc400/vsc40063/data/glen/gas-phase/H2
WORKDIR=$VSC_SCRATCH_NODE/$PBS_JOBID

echo Hostname: $(hostname)
echo ORIGDIR: $ORIGDIR
echo WORKDIR: $WORKDIR

mkdir -p $WORKDIR
cd $ORIGDIR
cp -r * $WORKDIR
ln -s $WORKDIR calc-dir
cd $WORKDIR
mympirun vasp > stdout
cp -r * $ORIGDIR
cd $ORIGDIR
rm calc-dir
rm -rf $WORKDIR
