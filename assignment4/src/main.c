/*
 * Copyright (C) NHR@FAU, University Erlangen-Nuremberg.
 * All rights reserved.
 * Use of this source code is governed by a MIT-style
 * license that can be found in the LICENSE file.
 */
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "likwid-marker.h"
#include "parameter.h"
#include "solver.h"
#include "timing.h"

int main(int argc, char** argv)
{
    double startTime, endTime;
    MPI_Init(&argc, &argv);
    Parameter params;
    Solver solver;
    initParameter(&params);

    if (argc < 2) {
        printf("Usage: %s <configFile>\n", argv[0]);
        exit(EXIT_SUCCESS);
    }
    readParameter(&params, argv[1]);
    printParameter(&params);

    initSolver(&solver, &params, 2);
    startTime = getTimeStamp();
    solveV3(&solver);
    endTime = getTimeStamp();

    getResult(&solver);

    if (solver.rank == 0) {
        printf("Walltime %.2fs\n", endTime - startTime);
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
}
