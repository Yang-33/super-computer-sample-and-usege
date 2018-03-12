#!/bin/bash
#PJM -L "rscgrp=fx-debug"
#PJM -L "node=12"
#PJM --mpi "proc=128"
#PJM -L "elapse=1:00"
fipp -C -d Prof mpirun ./wa2

# mpirun ./wa2



