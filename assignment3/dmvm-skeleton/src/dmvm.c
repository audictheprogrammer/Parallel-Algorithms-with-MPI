/*
 * Copyright (C) 2022 NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include "timing.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

double dmvm(double* restrict y,
    const double* restrict a,
    const double* restrict x,
    int N,
    int iter)
{
    double ts, te;

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rest = N % size;
    int N_local = N / size + (rank < rest);
    int lower_neighbor = (rank + 1) % size;
    int upper_neighbor = (rank - 1) + size*((rank-1) < 0);

    ts = MPI_Wtime();
    for (int j = 0; j < iter; j++) {

        for (int rot = 0; rot < size; rot++) {
            int cs = ((rank+rot) % size) * N_local; // column_shift, supposing N is equally distributed.
            for (int r = 0; r < N_local; r++) {
                for (int c = cs; c < cs+N_local; c++) {
                    y[r] = y[r] + a[r * N + c] * x[c-cs];
                }
            }

            MPI_Sendrecv_replace((void *)x, N_local, MPI_DOUBLE, upper_neighbor, 0, lower_neighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }
#ifdef CHECK
        {
            double sum = 0.0;

            for (int i = 0; i < N; i++) {
                sum += y[i];
                y[i] = 0.0;
            }
            fprintf(stderr, "Sum: %f\n", sum);
        }
#endif
    }
    te = MPI_Wtime();

    return te - ts;
}
