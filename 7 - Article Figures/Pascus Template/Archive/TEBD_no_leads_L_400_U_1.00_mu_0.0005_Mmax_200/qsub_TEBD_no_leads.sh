#!/bin/bash

cd /state/partition1

mkdir mocap_job$JOB_ID
cd mocap_job$JOB_ID


cp -r $SGE_O_WORKDIR/* ./

/share/apps/MATLAB/R2017a/bin/matlab  -nodisplay -nodesktop < ABTEBD_spinhalf_superfermion_no_leads_v1.m  >> out.dat

cp datafile.mat out.dat  $SGE_O_WORKDIR

cd ../

rm -R mocap_job$JOB_ID





