#!/bin/bash 

export OMP_NUM_THREADS=12

start_time=$(date +%s)
echo -e "\033[32m run ...\033[0m"

../../bin/MRCDFT -p 22Ne_para.dat -d 22Ne_b23.dat

echo calculation is finished !
end_time=$(date +%s)

execution_time=$((end_time - start_time))
execution_time_minutes=$((execution_time / 60))
execution_time_seconds=$((execution_time % 60))
echo "Time cost : ${execution_time_minutes}min${execution_time_seconds}s"

echo Done!



