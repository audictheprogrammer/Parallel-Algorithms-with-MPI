# Parallel Algorithms with MPI

This repository contains implementations and assignments focused on parallel programming using MPI (Message Passing Interface). The projects were completed as part of coursework at Friedrich-Alexander-Universität Erlangen-Nürnberg. Below is an overview of the repository and its content.

## Repository Structure

### Assignments

Each folder corresponds to an assignment. The descriptions and goals for each assignment are detailed below:

#### **Assignment 1: Warmup & Amdahl's Scaling Law**
- **Part 1**: Familiarization with the NHR@FAU cluster and development of a simple "Hello World" MPI application.
  - Implemented and tested a batch job for execution on Fritz nodes.
- **Part 2**: Analysis of Amdahl's Law for speedup measurements on multi-core systems.

#### **Assignment 2: Numerical Integration & Parallel Computation of π**
- **Part 2a**: Compute π using numerical integration with a rectangular integration scheme.
- **Part 2b**: Parallelize the numerical integration using MPI and implement a master/worker scheme for result aggregation.

#### **Assignment 3: Dense Matrix-Vector Multiplication (DMVM)**
- **Part (a)**: Perform strong scaling using blocking point-to-point communication.
- **Part (b)**: Optimize communication using non-blocking MPI calls and analyze performance improvements.

#### **Assignment 4: SOR Solver for the Poisson Problem**
- Parallel implementation of the Successive Over-Relaxation (SOR) method using a 1D domain decomposition.
- Conducted scalability studies and stability tests for varying domain sizes and process counts.
- Implemented a Red-Black SOR solver for better parallelization and compared it to the standard SOR solver.

#### **Assignment 5: Parallel CFD Solver with Derived Datatypes & Virtual Topologies**
- Parallelized a sequential CFD solver using:
  - Derived MPI datatypes.
  - 2D domain decomposition with virtual topologies.
- Conducted validation and benchmarking with test cases (lid-driven cavity and laminar channel flow).

#### **Assignment 6: Parallel 3D CFD Solver**
- Extended the concepts from Assignment 5 to implement a 3D parallel CFD solver with:
  - 3D domain decomposition.
  - Derived datatypes and MPI Cartesian topologies.
  - Neighborhood collective communication.
- Benchmarked performance on NFS and Lustre file systems using MPI I/O.

#### **Assignment 7: Benchmarking Pre-compiled Parallel 3D Solver**
- Conducted strong scaling studies for the provided parallelized 3D solver using:
  - Intra-NUMA domain scaling.
  - Inter-NUMA domain scaling.
  - Inter-node scaling.
- Analyzed communication, computation, and memory overhead based on scaling behavior.

## Tools and Technologies
- **Programming Languages**: C, MPI
- **HPC Resources**: Fritz cluster at NHR@FAU
- **Visualization**: Gnuplot, Paraview (for VTK file output)
- **Libraries and Utilities**: LIKWID for benchmarking, MPI I/O for file handling

## Getting Started

### Prerequisites
- Access to an HPC cluster with MPI installed.
- Modules used: Intel compiler, Intel MPI, and LIKWID.

### Compilation and Execution
- To compile the programs, use the provided `Makefile` in each assignment folder. For example:
  ```bash
  make
  ```
- To run the executables, use batch scripts or interactive job submissions. Example for submitting a job:
  ```bash
  sbatch job_script.sh
  ```

### Visualization
- Use Gnuplot for plotting results:
  ```bash
  gnuplot surface.plot
  ```
- View VTK files in Paraview:
  1. Open the file in Paraview.
  2. Click "Apply" to load the data.
  3. Use filters such as "Glyphs" to visualize velocity fields.

## Results and Observations
- Detailed results, including performance plots and conclusions, are documented in separate reports for each assignment.
- Performance metrics analyzed:
  - Wall-clock time.
  - Scalability (strong scaling).
  - Communication overhead.
  - Stability and convergence behavior.

## Acknowledgments
This work was completed as part of the "Parallel Algorithms" course at FAU. Special thanks to the instructors and the NHR@FAU HPC team for their support.

## License
This repository is private and intended for educational purposes. Unauthorized distribution or reproduction is prohibited.
