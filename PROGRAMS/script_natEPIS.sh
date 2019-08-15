#!/bin/bash
#$ -cwd

rm script_natEPIS.sh.*

#script_natEPIS.sh <LAMB> <NCRO> <BETA> <AVEs> <AVEh> <n>
#To be used only in a machine with /state/partition1 directory

#Check number of arguments
if [ $# -ne 6 ]  
then
	echo "Usage: $0 <LAMB> <NCRO> <BETA> <AVEs> <AVEh> <n>" 
	exit 1
fi

#Set arguments
LAMB=$1
NCRO=$2
BETA=$3
AVEs=$4
AVEh=$5
n=$6
 
#Working directory
WDIR=$PWD 

#Scratch directory
mkdir -p /state/partition1/sara$n/$SLURM_JOBID/

#Copy all files in scratch directory
cp seedfile /state/partition1/sara$n/$SLURM_JOBID/
cp ../naturalvEPIS /state/partition1/sara$n/$SLURM_JOBID/

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

#Move to scratch directory
cd /state/partition1/sara$n/$SLURM_JOBID

START=$(date +%s)
time ./naturalvEPIS>>out<<@
0
-99
1000	N
99	PS(99=random)
99	Lenght genome (99=free)
$NCRO	NCRO (max 2000)(Neu=Ncro)
30	NLOCI (2-30)
$LAMB	Lambda_a
0.015	Lambda_L
0.0	absolute effect of lethal (QT): normal (aL,aL)
0	random proportion(0) or large mutants (1)
1.0	Psi
$BETA	beta_s
$BETA	beta_a
$AVEs	ave |s|
$AVEs	ave |a|
0.0	PP_s
0.0	PP_a
2	dom model (0=cnt; 1:Deng, 2:CK94 gamma)
0.2	h_s (mod 0), k_s (mod 1)
$AVEh	ave h_s (mod 2)
0.2	h_a (mod 0), k_a (mod 1)
$AVEh	ave h_a (mod 2)
99	rho (99:a=s)
0	Vs
1	multi(1), add(2)
3	episK (1: no, 2: squared, 3: all)
2	episType (1: ncro, 2: whole)
0.0	sEPISmin
10000	generations
2000	gen/block
0	GENBOT
@

cat popfile >> popfile_L$LAMB.k$NCRO.s$AVEs.h$AVEh.B$BETA
cat datafile >> datafile_L$LAMB.k$NCRO.s$AVEs.h$AVEh.B$BETA
cat genfile >> genfile_L$LAMB.k$NCRO.s$AVEs.h$AVEh.B$BETA

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "naturalv took 		$DIFF seconds" >> timefile

#Copy output files to main directory
cp -r /state/partition1/sara$n/$SLURM_JOBID/popfile_L$LAMB.k$NCRO.s$AVEs.h$AVEh.B$BETA $WDIR/
cp -r /state/partition1/sara$n/$SLURM_JOBID/datafile_L$LAMB.k$NCRO.s$AVEs.h$AVEh.B$BETA $WDIR/
cp -r /state/partition1/sara$n/$SLURM_JOBID/genfile_L$LAMB.k$NCRO.s$AVEs.h$AVEh.B$BETA $WDIR/
cp -r /state/partition1/sara$n/$SLURM_JOBID/timefile $WDIR/

#Cleaning of scratch
rm -r /state/partition1/sara$n/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*
