#!/bin/bash

#SBATCH -N 1
#SBATCH --partition=genoa
#SBATCH --ntasks-per-node=192
#SBATCH --job-name=retracting_droplet_1st_sweep
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=a.s.bhargava@utwente.nl

MAXlevel="11"

Bo=( "0" "1.0" )

Oh=( "1E-2" "1E-1" "1.0" )

Gamma=( "4.0" "8.0" "16.0" )

theta=( "30.0" "90.0" "150.0" )

start=0

counter=0

for i in "${Bo[@]}"; do
    for j in "${Oh[@]}"; do

            dir=$(printf "%03d" $((start+counter)))
            mkdir -p $dir
 
            for k in "${Gamma[@]}"; do
                for l in "${theta[@]}"; do
                    cd $dir
                    srun -n 192 --gres=cpu:8 --exclusive retract_1st_sweep $MAXlevel $i $j $k $l & 
                    cd ../
                done
            done
            ((counter++))
            wait
    done
done

#CC99='mpicc -std=c99' qcc -Wall -O2 -disable-dimensions -D_MPI=1 retract_v2.c -o retract_1st_sweep -lm
