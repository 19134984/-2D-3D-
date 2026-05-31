#!/bin/bash

MPI_PROCS=3
OMP_THREADS_PER_NODE=24
OMP_THREADS_PER_RANK=$((OMP_THREADS_PER_NODE / MPI_PROCS))

export OMP_NUM_THREADS=${OMP_THREADS_PER_RANK}
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "======================================"
echo "Start time: $(date)"
echo "Host: $(hostname)"
echo "MPI processes: ${MPI_PROCS}"
echo "OpenMP threads per rank: ${OMP_NUM_THREADS}"
echo "Expected total threads: $((MPI_PROCS * OMP_NUM_THREADS))"
echo "======================================"

mpirun -np ${MPI_PROCS} ./a.exe

echo "======================================"
echo "End time: $(date)"
echo "======================================"