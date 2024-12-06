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

#include "affinity.h"
#include "allocate.h"
#include "timing.h"

int main(int argc, char** argv) {
    if (argc == 1) {
        printf("Buongiorno Mondo \n");
        return EXIT_SUCCESS; 
    }
    printf("Buongiorno %s \n", argv[1]);
    return EXIT_SUCCESS; 
}
