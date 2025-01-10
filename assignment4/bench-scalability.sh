#!/bin/bash -l
#SBATCH --job-name=bench_scalability_seq_24_1000
#SBATCH --output=scalability_results_%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=6:00:00
#SBATCH --cpu-freq=2400000-2400000:performance
#SBATCH --exclusive

module load intel intelmpi

RESULTS_FILE="scalability_results_seq_omega1.2_cpu2400000_SIZE1000.csv"
RAW_OUTPUT_DIR="raw_outputs"

mkdir -p $RAW_OUTPUT_DIR

echo "tasks,iterations,walltime,status" > $RESULTS_FILE

# Fixed parameters.
XLENGTH=1.0
YLENGTH=1.0
IMAX=1000
JMAX=1000
ITERMAX=10000000
EPS=0.000001

for NB_TASK in $(seq 1 32); do
    echo "Running scalability test: nb_tasks=$NB_TASK"
    RAW_OUTPUT_FILE="$RAW_OUTPUT_DIR/output_n_$NB_TASK.txt"

    # Running the task
    mpirun -np $NB_TASK ./exe-ICC poisson.par > $RAW_OUTPUT_FILE 2>&1

    # Extract values directly from the raw output file.
    ITERATIONS=$(grep "took" $RAW_OUTPUT_FILE | awk '{print $5}')
    WALLTIME=$(grep "Walltime" $RAW_OUTPUT_FILE | awk '{print $2}')
    STATUS=$(grep "took" $RAW_OUTPUT_FILE | awk '{print $8}')

    # Default values if parsing fails.
    ITERATIONS=${ITERATIONS:-"N/A"}
    WALLTIME=${WALLTIME:-"N/A"}
    STATUS=${STATUS:-"N/A"}

    # Append the result to the CSV file.
    echo "$NB_TASK,$ITERATIONS,$WALLTIME,$STATUS" >> $RESULTS_FILE
done

echo "Scalability analysis completed. Results saved to $RESULTS_FILE and raw outputs to $RAW_OUTPUT_DIR."
