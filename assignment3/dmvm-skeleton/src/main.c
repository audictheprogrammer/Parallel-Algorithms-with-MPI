/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>

#include "allocate.h"
#include "timing.h"

extern double dmvm(double* restrict y,
    const double* restrict a,
    const double* restrict x,
    int N,
    int iter);

int main(int argc, char** argv)
{

    size_t bytesPerWord = sizeof(double);
    size_t N            = 0;
    size_t iter         = 1;
    double *a, *x, *y;
    double t0, t1;
    double walltime;

    if (argc > 2) {
        N    = atoi(argv[1]);
        iter = atoi(argv[2]);
    } else {
        printf("Usage: %s <N> <iter>\n", argv[0]);
        exit(EXIT_SUCCESS);
    }

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rest = N % size;
    int N_local = N / size + (rank < rest);
    int N_allocate = N / size + rest;

    a = (double*)allocate(ARRAY_ALIGNMENT, N_local * N * bytesPerWord);
    x = (double*)allocate(ARRAY_ALIGNMENT, N_allocate * bytesPerWord);
    y = (double*)allocate(ARRAY_ALIGNMENT, N_local * bytesPerWord);

    // initialize arrays
    for (int i = 0; i < N_local; i++) {
        double I = rank * (N/size) + fmin(rank, rest) + i; // Global I vs Local i. Assuming N % size == 0.
        x[i] = (double)I;
        y[i] = 0.0;

        for (int j = 0; j < N; j++) {
            a[i * N + j] = (double)j + I;
        }
    }

    walltime = dmvm(y, a, x, N, iter);
    free(x);
    free(y);
    free(a);
    double flops = (double)2.0 * N * N * iter;
    // # iterations, problem size, flop rate, walltime
    if (rank == 0) {
        printf("\n%zu %zu %.2f %.2f\n", iter, N, 1.0E-06 * flops / walltime, walltime);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
