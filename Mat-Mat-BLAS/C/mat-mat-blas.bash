#!/bin/bash
#PJM -L "rscgrp=fx-debug"
#PJM -L "node=1"
#PJM --mpi "proc=1"
#PJM -L "elapse=1:00"
export OMP_NUM_THREADS=32
mpirun ./mat-mat-blas



