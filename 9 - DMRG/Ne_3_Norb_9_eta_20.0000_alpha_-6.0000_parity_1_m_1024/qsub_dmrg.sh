#!/bin/bash


dir_sf=$SGE_O_WORKDIR


cd /state/partition1

mkdir wigner_job$JOB_ID
cd wigner_job$JOB_ID


cp $SGE_O_WORKDIR/* ./

/share/apps/MATLAB/R2017a/bin/matlab  -nodisplay -nodesktop < wc_project.m  >> out.dat


#cp FCI* $$SGE_O_WORKDIR/
#cp classical_positions_*.dat $$SGE_O_WORKDIR/
#cd N4


cp -r  *  $SGE_O_WORKDIR/

cd ../

rm -R wigner_job$JOB_ID




