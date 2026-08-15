#!/usr/bin/env bash
set -e
export MODULEPATH=/usr/share/Modules/modulefiles
source /scratch/apps/Modules/init/bash
module purge
module load mpi/latest
module load fftw/3.3.10-intel-2023
module load hdf5/1.10.6-intel-2023
source /share/apps/intel/oneapi/2023.2.0/setvars.sh
cd "/share/org/SHUTUANL/shu_zhangjs/rainy model/rotating_case/latest_program/source"
make obj/force_balance_online.o obj/openfi.o obj/gcurv.o simexec
