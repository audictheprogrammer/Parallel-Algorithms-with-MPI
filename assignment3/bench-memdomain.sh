#!/bin/bash -l
#SBATCH --job-name=bench_memdomain
#SBATCH --output=Fritz_ICX_DMVM_memdomain
#SBATCH --partition=singlenode
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --time=10:00:00
#SBATCH --export=NONE
#SBATCH --cpu-freq=2400000-2400000:performance

unset SLURM_EXPORT_ENV

module load likwid intel intelmpi
export I_MPI_PIN=1
export I_MPI_DEBUG=0

FILENAME="result_bench_memdomain.csv"

cd ~/Parallel-Algorithms-with-MPI/assignment3/dmvm-skeleton/
make distclean
make

rm -f $FILENAME
touch $FILENAME
echo "Ranks,NITER,N,MFlops,Time" >>$FILENAME

_iterate() {
    for np in $(seq 1 $NPM); do
        np_1=$(($np - 1))
        export I_MPI_PIN_PROCESSOR_LIST=0-$np_1

        raw_result=$(mpirun -n $np ./exe-ICX $N $NITER 2>&1)
        echo "Running mpirun -n $np ./exe-ICX $N $NITER"
        
        result=$(echo "$raw_result" | grep -E '^[0-9]+ [0-9]+ [0-9.]+ [0-9.]+$')

        if [[ ! -z $result ]]; then
            echo "$np $result" >>$FILENAME
        else
            echo "Warning: No valid result captured for np=$np, NITER=$NITER, N=$N" >> debug.log
        fi
    done
}

NPM=18

# For domain of 1000x1000
NITER=1000000
N=1000
_iterate

# For domain of 4000x4000
NITER=100000
N=4000
_iterate

# For domain of 10000x10000
NITER=10000
N=10000
_iterate

# For domain of 20000x20000
NITER=5000
N=20000
_iterate

sed -i 's/ /,/g' $FILENAME
