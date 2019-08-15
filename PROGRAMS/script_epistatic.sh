#!/bin/bash
#$ -cwd

rm script_epistatic.sh.*

#To be used only in a machine with /scratch directory

#Check number of arguments
if [ $# -ne 1 ]  
then
	echo "Usage: $0 <n>" 
	exit 1
fi

#Set arguments
n=$1

#Working directory
WDIR=$PWD 

#Scratch directory
mkdir -p /state/partition1/epistaticEepis$n/$SLURM_JOBID/

#Copy all files in scratch directory
cp seedfile /state/partition1/epistaticEepis$n/$SLURM_JOBID/
cp datafile* /state/partition1/epistaticEepis$n/$SLURM_JOBID/
cp popfile* /state/partition1/epistaticEepis$n/$SLURM_JOBID/
cp ../epistatic /state/partition1/epistaticEepis$n/$SLURM_JOBID/

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

#Move to scratch directory
cd /state/partition1/epistaticEepis$n/$SLURM_JOBID

########### CASOS ########### 

NCRO=100
REP=1000
LAMB=0.05

for AVEs in 0.05 0.1 0.2
do
for AVEh in 0.2 0.3 0.4
do
for BETA in 2.0
do

rm popfile
rm datafile

cp popfile_L$LAMB.k$NCRO.s$AVEs.h$AVEh.B$BETA popfile 
cp datafile_L$LAMB.k$NCRO.s$AVEs.h$AVEh.B$BETA datafile 


#Run the programme epistatic.c type=0 PS=99


START=$(date +%s)
time ./epistatic>>out<<@
0
-99
2	NMAX
42	MAXFEC
11	DS_MAXFEC
99	PS (Probability of selfing; 99=random)
99	L
$NCRO	NCRO (max 2000)(Neu=Ncro)
30	NLOCI (2-30)
$LAMB	Lambda_s
0.015	Lambda_L
$BETA	beta_s
$AVEs	ave |s|
1	dom (0:constant; 1:variable)
$AVEh	ave h_s
0	relaxation factor (0: no, 1:yes)
0	neutral (0: no, 1:yes)
1	scaling (0: no, 1:yes)
3	episK (1: no, 2: squared, 3: all)
2	episType (1: ncro, 2: whole)
0.0	sEPISmin
7	generations
$REP	replicates
@
cat genfile.dat >> GENFILEext
cat salida.dat >> SALIDA
#cat outfile.dat >> OUTFILE

cp -r /state/partition1/epistaticEepis$n/$SLURM_JOBID/GENFILEext $WDIR/
cp -r /state/partition1/epistaticEepis$n/$SLURM_JOBID/SALIDA $WDIR/
cp -r /state/partition1/epistaticEepis$n/$SLURM_JOBID/timefile $WDIR/

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "extinction took 	$DIFF seconds" >> timefile

done
done
done

#Copy output files to main directory
cp -r /state/partition1/epistaticEepis$n/$SLURM_JOBID/seedfile $WDIR
cp -r /state/partition1/epistaticEepis$n/$SLURM_JOBID/GENFILEext $WDIR/
cp -r /state/partition1/epistaticEepis$n/$SLURM_JOBID/timefile $WDIR/
#cp -r /state/partition1/epistaticEepis$n/$SLURM_JOBID/pedfile.dat $WDIR/
cp -r /state/partition1/epistaticEepis$n/$SLURM_JOBID/dfilename.dat $WDIR/
#cp -r /state/partition1/epistaticEepis$n/$SLURM_JOBID/Wreplicates $WDIR/
cp -r /state/partition1/epistaticEepis$n/$SLURM_JOBID/SALIDA $WDIR/

#Cleaning of scratch
rm -r /state/partition1/epistaticEepis$n/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*
