#!/bin/bash 
#PBS -o 150Nd.out
#PBS -e 150Nd.err
#PBS -N 150Nd
#PBS -l select=1:ncpus=64:host=cn5

# export OMP_NUM_THREADS=8
# export MKL_NUM_THREADS=8
# export MKL_DYNAMIC=FALSE
# export MKL_THREADING_LAYER=GNU

start_time=$(date +%s)
echo -e "\033[32m run ...\033[0m"

cd /home/yli/MR_CDFT/MRCDFT_dev/MRCDFT_sens/examples/150Nd
mpiexec -np 10 ../../bin/MRCDFT -p 150Nd_para.dat -d 150Nd_b23.dat -f 150Nd_forcefactors.dat

echo calculation is finished !
end_time=$(date +%s)
execution_time=$((end_time - start_time))
execution_time_minutes=$((execution_time / 60))
execution_time_seconds=$((execution_time % 60))
echo "Time cost : ${execution_time_minutes}min${execution_time_seconds}s"
echo Done!



