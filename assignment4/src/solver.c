/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved. This file is part of nusif-solver.
 * Use of this source code is governed by a MIT style
 * license that can be found in the LICENSE file.
 */
#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include "allocate.h"
#include "parameter.h"
#include "solver.h"

#define PI        3.14159265358979323846
#define P(i, j)   p[(j) * (imax + 2) + (i)]
#define RHS(i, j) rhs[(j) * (imax + 2) + (i)]

static int sizeOfRank(int rank, int size, int N)
{
    return (N / size) + (rank < (N % size));
}

static void exchange(Solver* solver)
{
    MPI_Request requests[4] = { MPI_REQUEST_NULL,
        MPI_REQUEST_NULL,
        MPI_REQUEST_NULL,
        MPI_REQUEST_NULL };

    /* exchange ghost cells with top neighbor */
    if (solver->rank + 1 < solver->size) {
        int top     = solver->rank + 1;
        double* src = &(solver->p[(solver->imax + 2) * (solver->jmaxLocal + 0)]);
        double* dst = &(solver->p[(solver->imax + 2) * (solver->jmaxLocal + 1)]);

        MPI_Isend(src, solver->imax + 2, MPI_DOUBLE, top, 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(dst, solver->imax + 2, MPI_DOUBLE, top, 1, MPI_COMM_WORLD, &requests[1]);
    }

    /* exchange ghost cells with bottom neighbor */
    if (solver->rank > 0) {
        int bottom  = solver->rank - 1;
        double* src = &(solver->p[solver->imax + 2]);
        double* dst = &(solver->p[0]);

        MPI_Isend(src, solver->imax + 2, MPI_DOUBLE, bottom, 1, MPI_COMM_WORLD, &requests[2]);
        MPI_Irecv(dst, solver->imax + 2, MPI_DOUBLE, bottom, 0, MPI_COMM_WORLD, &requests[3]);
    }

    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
}

void getResult(Solver* solver)
{
    if (solver->size == 1) {
        writeResult(solver, "p.dat");
        return ;
    }
    double* p = NULL;
    int *rcvCounts, *displs;

    if (solver->rank == 0) {
        p = allocate(64, (solver->imax + 2) * (solver->jmax + 2) * sizeof(double));
        rcvCounts    = (int*) malloc(solver->size * sizeof(int));
        displs       = (int*) malloc(solver->size * sizeof(int));
        rcvCounts[0] = (solver->jmaxLocal + 1) * (solver->imax + 2);
        displs[0]    = 0;
        int cursor   = rcvCounts[0];

        for (int i = 1; i < solver->size - 1; i++) {
            rcvCounts[i] = solver->jmaxLocal * (solver->imax + 2);
            displs[i]    = cursor;
            cursor += rcvCounts[i];
        }

        rcvCounts[solver->size - 1] = (solver->jmaxLocal + 1) * (solver->imax + 2);
        displs[solver->size - 1] = cursor;
    }

    int cnt = (solver->jmaxLocal) * (solver->imax + 2);
    if (solver->rank == 0 || solver->rank == solver->size - 1) {
        cnt += solver->imax + 2;
    }

    double* sendbuffer = &(solver->p[solver->imax + 2]);
    if (solver->rank == 0) {
        sendbuffer = &(solver->p[0]);
    }

    MPI_Gatherv(sendbuffer, cnt, MPI_DOUBLE, p, rcvCounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Memory free.
    // free(solver->p);
    // free(solver->rhs);

    if (solver->rank == 0) {
        solver->p = p;
        solver->jmaxLocal = solver->jmax;
        writeResult(solver, "p.dat");

        // Memory free
        // free(displs);
        // free(rcvCounts);
        // free(p);
    }
}

void initSolver(Solver* solver, Parameter* params, int problem)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &(solver->rank));
    MPI_Comm_size(MPI_COMM_WORLD, &(solver->size));
    solver->imax = params->imax;
    solver->jmax = params->jmax;
    solver->jmaxLocal = sizeOfRank(solver->rank, solver->size, solver->jmax);
    printf("RANK %d: imaxLocal : %d, jmaxLocal : %d\n",
        solver->rank,
        solver->imax,
        solver->jmaxLocal);

    solver->dx = params->xlength / params->imax;
    solver->dy = params->ylength / params->jmax;
    // solver->ys      = solver->rank * solver->jmaxLocal * solver->dy;
    solver->eps     = params->eps;
    solver->omega   = params->omg;
    solver->itermax = params->itermax;

    int rank      = solver->rank;
    int size      = solver->size;
    int imax      = solver->imax;
    int jmax      = solver->jmax;
    int jmaxLocal = solver->jmaxLocal;
    solver->p     = allocate(64, (imax + 2) * (jmaxLocal + 2) * sizeof(double));
    solver->rhs   = allocate(64, (imax + 2) * (jmax + 2) * sizeof(double));

    double dx   = solver->dx;
    double dy   = solver->dy;
    double* p   = solver->p;
    double* rhs = solver->rhs;

    for (int j = 0; j < jmaxLocal + 2; j++) {
        int J = rank * (jmax / size) + MIN(rank, jmax % size) + j;
        for (int i = 0; i < imax + 2; i++) {
            P(i, j) = sin(2.0 * PI * i * dx * 2.0) + sin(2.0 * PI * J * dy * 2.0);
        }
    }

    if (problem == 2) {
        for (int j = 0; j < jmax + 2; j++) {
            for (int i = 0; i < imax + 2; i++) {
                RHS(i, j) = sin(2.0 * PI * i * dx);
            }
        }
    } else {
        for (int j = 0; j < jmax + 2; j++) {
            for (int i = 0; i < imax + 2; i++) {
                RHS(i, j) = 0.0;
            }
        }
    }
}

/* Solve V1: u^{t+1}_{i, j} = F(u^{t+1}_{i-1, j},
                                   u^{t+1}_{i, j-1},
                                   u^{t}_{i+1, j},
                                   u^{t}_{i, j+1}).
*/
void solve(Solver* solver)
{
    int imax      = solver->imax;
    int jmax      = solver->jmax;
    int jmaxLocal = solver->jmaxLocal;
    double eps    = solver->eps;
    int itermax   = solver->itermax;
    double dx2    = solver->dx * solver->dx;
    double dy2    = solver->dy * solver->dy;
    double idx2   = 1.0 / dx2;
    double idy2   = 1.0 / dy2;
    double factor = solver->omega * 0.5 * (dx2 * dy2) / (dx2 + dy2);
    double* p     = solver->p;
    double* rhs   = solver->rhs;
    double epssq  = eps * eps;
    int it        = 0;
    double res    = eps + 1.0;

    while ((res >= epssq) && (it < itermax)) {
        res = 0.0;
        exchange(solver);

        for (int j = 1; j < jmaxLocal + 1; j++) {
            for (int i = 1; i < imax + 1; i++) {
                double r = RHS(i, j) -
                           ((P(i - 1, j) - 2.0 * P(i, j) + P(i + 1, j)) * idx2 +
                               (P(i, j - 1) - 2.0 * P(i, j) + P(i, j + 1)) * idy2);

                P(i, j) -= (factor * r);
                res += (r * r);
            }
        }

        if (solver->rank == 0) {
            for (int i = 1; i < imax + 1; i++) {
                P(i, 0) = P(i, 1);
            }
        }
        if (solver->rank == solver->size - 1) {
            for (int i = 1; i < imax + 1; i++) {
                P(i, jmaxLocal + 1) = P(i, jmaxLocal);
            }
        }

        for (int j = 1; j < jmaxLocal + 1; j++) {
            P(0, j)        = P(1, j);
            P(imax + 1, j) = P(imax, j);
        }

        MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        res = res / (double)(imax * jmax);
#ifdef DEBUG
        // if (solver->rank == 0) {
        printf("%d Residuum: %e\n", it, res);
        // }
#endif
        it++;
    }

    if (solver->rank == 0) {
        printf("Solver Gauss-Seidel: took %d iterations to reach %f\n",
                it,
                sqrt(res));
        }
}

/* Solve V2 Red-Black Gauss-Seidel: u^{RED}_{i, j} = F(u^{BLACK}_{i-1, j},
                                   u^{BLACK}_{i, j-1},
                                   u^{BLACK}_{i+1, j},
                                   u^{BLACK}_{i, j+1}).
Same for u^{BLACK}_{i, j}.
*/
void solveV2(Solver* solver)
{
    int imax      = solver->imax;
    int jmax      = solver->jmax;
    int jmaxLocal = solver->jmaxLocal;
    int rank      = solver->rank;
    int size      = solver->size;
    double eps    = solver->eps;
    int itermax   = solver->itermax;
    double dx2    = solver->dx * solver->dx;
    double dy2    = solver->dy * solver->dy;
    double idx2   = 1.0 / dx2;
    double idy2   = 1.0 / dy2;
    double factor = solver->omega * 0.5 * (dx2 * dy2) / (dx2 + dy2);
    double* p     = solver->p;
    double* rhs   = solver->rhs;
    double epssq  = eps * eps;
    int it        = 0;
    double res    = eps + 1.0;

    while ((res >= epssq) && (it < itermax)) {
        res = 0.0;
        exchange(solver);

        // Red cells only.
        for (int j = 1; j < jmaxLocal + 1; j++) {
            int J = rank * (jmax / size) + MIN(rank, jmax % size) + j;
            for (int i = 1; i < imax + 1; i++) {
                if ((i+J) % 2 == 0) {
                    continue;
                }
                double r = RHS(i, j) -
                           ((P(i - 1, j) - 2.0 * P(i, j) + P(i + 1, j)) * idx2 +
                               (P(i, j - 1) - 2.0 * P(i, j) + P(i, j + 1)) * idy2);

                P(i, j) -= (factor * r);
                res += (r * r);
            }
        }

        exchange(solver);

        // Black cells only.
        for (int j = 1; j < jmaxLocal + 1; j++) {
            int J = rank * (jmax / size) + MIN(rank, jmax % size) + j;
            for (int i = 1; i < imax + 1; i++) {
                if ((i+J) % 2 == 1) {
                    continue;
                }
                double r = RHS(i, j) -
                           ((P(i - 1, j) - 2.0 * P(i, j) + P(i + 1, j)) * idx2 +
                               (P(i, j - 1) - 2.0 * P(i, j) + P(i, j + 1)) * idy2);

                P(i, j) -= (factor * r);
                res += (r * r);
            }
        }

        if (solver->rank == 0) {
            for (int i = 1; i < imax + 1; i++) {
                P(i, 0) = P(i, 1);
            }
        }
        if (solver->rank == solver->size - 1) {
            for (int i = 1; i < imax + 1; i++) {
                P(i, jmaxLocal + 1) = P(i, jmaxLocal);
            }
        }

        for (int j = 1; j < jmaxLocal + 1; j++) {
            P(0, j)        = P(1, j);
            P(imax + 1, j) = P(imax, j);
        }

        MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        res = res / (double)(imax * jmax);
#ifdef DEBUG
        // if (solver->rank == 0) {
        printf("%d Residuum: %e\n", it, res);
        // }
#endif
        it++;
    }

    if (solver->rank == 0) {
        printf("Solver Gauss-Seidel REDBLACK: took %d iterations to reach %f\n",
                it,
                sqrt(res));
        }
}


/* Solve V2 Red-Black SOR: u^{RED}_{i, j} = F(u^{BLACK}_{i-1, j},
                                   u^{BLACK}_{i, j-1},
                                   u^{BLACK}_{i+1, j},
                                   u^{BLACK}_{i, j+1},
                                   omega).
Same for u^{BLACK}_{i, j}.
*/
void solveV3(Solver* solver)
{
    int imax      = solver->imax;
    int jmax      = solver->jmax;
    int jmaxLocal = solver->jmaxLocal;
    int rank      = solver->rank;
    int size      = solver->size;
    double eps    = solver->eps;
    int itermax   = solver->itermax;
    double dx2    = solver->dx * solver->dx;
    double dy2    = solver->dy * solver->dy;
    double idx2   = 1.0 / dx2;
    double idy2   = 1.0 / dy2;
    double omega  = solver->omega;
    double factor = omega * 0.5 * (dx2 * dy2) / (dx2 + dy2);
    double* p     = solver->p;
    double* rhs   = solver->rhs;
    double epssq  = eps * eps;
    int it        = 0;
    double res    = eps + 1.0;

    while ((res >= epssq) && (it < itermax)) {
        res = 0.0;
        exchange(solver);

        // Red cells only.
        for (int j = 1; j < jmaxLocal + 1; j++) {
            int J = rank * (jmax / size) + MIN(rank, jmax % size) + j;
            for (int i = 1; i < imax + 1; i++) {
                if ((i+J) % 2 == 0) {
                    continue;
                }
                double r = RHS(i, j) -
                           ((P(i - 1, j) - 2.0 * P(i, j) + P(i + 1, j)) * idx2 +
                               (P(i, j - 1) - 2.0 * P(i, j) + P(i, j + 1)) * idy2);

                P(i, j) -= (omega * factor * r);
                res += (r * r);
            }
        }

        exchange(solver);

        // Black cells only.
        for (int j = 1; j < jmaxLocal + 1; j++) {
            int J = rank * (jmax / size) + MIN(rank, jmax % size) + j;
            for (int i = 1; i < imax + 1; i++) {
                if ((i+J) % 2 == 1) {
                    continue;
                }
                double r = RHS(i, j) -
                           ((P(i - 1, j) - 2.0 * P(i, j) + P(i + 1, j)) * idx2 +
                               (P(i, j - 1) - 2.0 * P(i, j) + P(i, j + 1)) * idy2);

                P(i, j) -= (omega * factor * r);
                res += (r * r);
            }
        }

        if (solver->rank == 0) {
            for (int i = 1; i < imax + 1; i++) {
                P(i, 0) = P(i, 1);
            }
        }
        if (solver->rank == solver->size - 1) {
            for (int i = 1; i < imax + 1; i++) {
                P(i, jmaxLocal + 1) = P(i, jmaxLocal);
            }
        }

        for (int j = 1; j < jmaxLocal + 1; j++) {
            P(0, j)        = P(1, j);
            P(imax + 1, j) = P(imax, j);
        }

        MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        res = res / (double)(imax * jmax);
#ifdef DEBUG
        // if (solver->rank == 0) {
        printf("%d Residuum: %e\n", it, res);
        // }
#endif
        it++;
    }

    if (solver->rank == 0) {
        if (isnan(res)) {
            printf("Solver SOR REDBLACK: took %d iterations to DV using omega=%f\n",
                    it,
                    solver->omega);
        } else {
            printf("Solver SOR REDBLACK: took %d iterations to CV using omega=%f\n",
                    it,
                    solver->omega);
        }
    }
}

void writeResult(Solver* solver, char* filename)
{
    int imax      = solver->imax;
    int jmax      = solver->jmax;
    int jmaxLocal = solver->jmaxLocal;
    double* p     = solver->p;

    FILE* fp;
    fp = fopen(filename, "w");

    if (fp == NULL) {
        printf("Error!\n");
        exit(EXIT_FAILURE);
    }

    for (int j = 0; j < jmaxLocal + 2; j++) {
        for (int i = 0; i < imax + 2; i++) {
            fprintf(fp, "%f ", P(i, j));
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}
