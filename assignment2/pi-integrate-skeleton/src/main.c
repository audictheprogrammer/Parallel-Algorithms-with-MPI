/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <mpi.h>

#include "timing.h"

double integrate(double, double);

int main (int argc, char** argv) {
    double wcs, wce;
    double Pi;
    double a = 0, b = 1;
    
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double a_local = a + rank * ((b-a) / size);
    double b_local = a + (rank+1) * ((b-a) / size);

    wcs = MPI_Wtime();
    double Pi_loc = integrate(a_local, b_local);
    MPI_Allreduce(&Pi_loc, &Pi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    Pi *= 4;
    wce = MPI_Wtime();

    printf("Pi=%.15lf in %.3lf s \n", Pi,wce-wcs);
    MPI_Finalize();
    return EXIT_SUCCESS;
}

double f(double x) {
    return sqrt(1 - pow(x, 2));
}

double integrate(double a, double b) {
	/* Your logic to integrate between given interval a to b.
    Declare N here and calculate h using a, b and N.
    Iterate over N, calculate the area and sum them.
    Return sum * h to get the value of PI. */
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = 1e9;
    double N_local = N / size;
    double sum = 0.;
    double h = (b-a) / N_local;

    for (int i = 0; i < N_local; i++) {
        sum += f(a + i*h);
    }

	return h*sum;
}
