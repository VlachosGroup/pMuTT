#!/bin/bash
#$ -cwd
#$ -pe openmpi-smp 16
#$ -q *@@3rd_gen
#$ -j y
#$ -S /bin/bash
#
# Get our environment setup:
#
source /etc/profile.d/valet.sh
vpkg_require "vasp/5.3.2+vtst+d3-pgi10"
#
INPUTFILES="POSCAR POTCAR INCAR KPOINTS"
#
OUTPUTFILES="CONTCAR OUTCAR OSZICAR"
#
## You should NOT need to modify anything after this comment!!
#
Cleanup()
{
  for ofile in $OUTPUTFILES; do
    if [ -e $TMPDIR/${ofile} ]; then
      cp -p $TMPDIR/${ofile} $SGE_O_WORKDIR/${ofile}
    fi
  done
 
  if [ $# -eq 1 ] ; then
    exit $1
  fi
  exit 0
}
# Copy the input files into place:
for ifile in $INPUTFILES; do
  if [ -e $SGE_O_WORKDIR/${ifile} ]; then
    cp -rf $SGE_O_WORKDIR/${ifile} $TMPDIR/${ofile}
  fi
done
echo "copy finished"
 
# Sync the temp directory to the other nodes:
NODES=`cat $PE_HOSTFILE|cut -f1 -d" "|cut -f1 -d"."`
for node in $NODES; do
  if [ ${HOSTNAME} != ${node} ]; then
    echo "Replicating $TMPDIR to $node"
    rcp -r $TMPDIR ${node}:$TMPDIR
  fi
done
echo "sync finished"
 
# Begin the run by printing some basic info and then
# invoke mpiexec:
echo "GridEngine parameters:"
echo "  nhosts         = $NHOSTS"
echo "  nproc          = $NSLOTS"
echo "  mpiexec        =" `which mpiexec`
echo "  vasp           =" `which vasp`
echo "-- begin OPENMPI run --"
cd $TMPDIR
time mpiexec --mca btl sm,self --display-map vasp
Cleanup
echo "-- end OPENMPI run --"
