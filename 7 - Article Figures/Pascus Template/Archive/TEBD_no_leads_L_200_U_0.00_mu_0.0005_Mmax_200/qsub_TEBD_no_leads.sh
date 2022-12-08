#!/bin/bash

cd /state/partition1

mkdir mocap_job$JOB_ID
cd mocap_job$JOB_ID

cp -r $SGE_O_WORKDIR/* ./

  echo Setting up environment variables
  MCRROOT="/share/apps/MATLAB/R2022b"
  echo ---
  LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
  export LD_LIBRARY_PATH;
  echo LD_LIBRARY_PATH is ${LD_LIBRARY_PATH};
# Preload glibc_shim in case of RHEL7 variants
  test -e /usr/bin/ldd &&  ldd --version |  grep -q "(GNU libc) 2\.17"  \
            && export LD_PRELOAD="${MCRROOT}/bin/glnxa64/glibc-2.17_shim.so"
  shift 1
 

exe_name=ABTEBD_spinhalf_superfermion_no_leads_v2
eval "\"./ABTEBD_spinhalf_superfermion_no_leads_v2\""  > out.txt

cp datafile.mat out.dat  $SGE_O_WORKDIR

cd ../

rm -R mocap_job$JOB_ID





