#!/bin/bash
#$ -cwd
#$ -N CH4
#$ -pe openmpi-smp 16
#$ -q *@@3rd_gen
#$ -j y
#$ -S /bin/bash
#
# Get our environment setup:
#
source /etc/profile.d/valet.sh
vpkg_require "vasp/5.3.2+vtst+d3+gamma-pgi10"
#
## You should NOT need to modify anything after this comment!!
#
 
# Begin the run by printing some basic info and then
# invoke mpiexec:
echo "GridEngine parameters:"	
echo "  nhosts         = $NHOSTS" 
echo "  nproc          = $NSLOTS" 
echo "  mpiexec        =" `which mpiexec` 
echo "  pe_hostfile    = $PE_HOSTFILE" 
echo "  vasp           =" `which vasp`	
echo 
cat $PE_HOSTFILE 
echo 
echo "-- begin OPENMPI run --"
time mpiexec --n $NSLOTS --host localhost --mca btl sm,self --display-map vasp
echo "-- end OPENMPI run --"
