#!/bin/bash -l
#SBATCH --job-name=bench_convergence_seq_1000
#SBATCH --output=convergence_results_%j.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --cpu-freq=2400000-2400000:performance
#SBATCH --exclusive

module load intel intelmpi

# Fixed parameters.
XLENGTH=1.0
YLENGTH=1.0
IMAX=1000
JMAX=1000
ITERMAX=10000000
EPS=0.000001

RESULTS_FILE="convergence_seq_results_$JMAX.csv"
RAW_OUTPUT_DIR="raw_outputs"

mkdir -p $RAW_OUTPUT_DIR

echo "omega,iterations,walltime,status" > $RESULTS_FILE


# Omega values.
OMEGAS=($(seq 1.00 0.1 1.99))

for OMG in "${OMEGAS[@]}"; do
    # Generating a temporary poisson.par file.
    TEMP_PAR="poisson_$OMG.par"
    cat <<EOL > $TEMP_PAR
name poisson

xlength    $XLENGTH
ylength    $YLENGTH
imax       $IMAX
jmax       $JMAX

itermax  $ITERMAX
eps      $EPS
omg      $OMG
EOL

    echo "Running with omega=$OMG"
    RAW_OUTPUT_FILE="$RAW_OUTPUT_DIR/output_omega_$OMG.txt"

    mpirun -n 1 ./exe-ICC $TEMP_PAR > $RAW_OUTPUT_FILE 2>&1

    # Extract values directly from the raw output file.
    ITERATIONS=$(grep "took" $RAW_OUTPUT_FILE | awk '{print $5}')
    WALLTIME=$(grep "Walltime" $RAW_OUTPUT_FILE | awk '{print $2}')
    STATUS=$(grep "took" $RAW_OUTPUT_FILE | awk '{print $8}')

    # Default values if parsing fails.
    ITERATIONS=${ITERATIONS:-"N/A"}
    WALLTIME=${WALLTIME:-"N/A"}
    STATUS=${STATUS:-"N/A"}

    echo "$OMG,$ITERATIONS,$WALLTIME,$STATUS" >> $RESULTS_FILE

    rm -f $TEMP_PAR
done

echo "Convergence analysis completed. Results saved to $RESULTS_FILE and raw outputs to $RAW_OUTPUT_DIR."
