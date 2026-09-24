#!/bin/bash

# Recommended for OpenMP tasks

#SBATCH -J YW-job             # Job name
#SBATCH -o job.%j.out         # Name of stdout output file (%j expands to jobId)
#SBATCH -e stderr             # Name of stderr output file (%j expands to jobId)
#SBATCH -N 1                  # Request that a minimum of minnodes nodes be allocated to this job. A maximum node count may also be specified with maxnodes. If only one number is specified, this is used as both the minimum and maximum node count. If -N is not specified, the default behavior is to allocate enough nodes to satisfy the requested resources as expressed by per-job specification options, e.g. -n, -c and --gpus.
#SBATCH -n 1                  # Total number of MPI tasks (processes)
#SBATCH -c 64                 # cpus-per-task (process) (will be placed in one node, but not necessarily in one socket)
#SBATCH --sockets-per-node=1   # Restrict node selection to nodes with at least the specified number of sockets. (only as a lower bound!)
#SBATCH --cores-per-socket=64  # Restrict node selection to nodes with at least the specified number of cores per socket. (only as a lower bound!)
#SBATCH --threads-per-core=1   # Restrict node selection to nodes with at least the specified number of threads per core. In task layout, use the specified maximum number of threads per core. (exact value!)
#SBATCH --mem=16G             # memory per node
#SBATCH -t 16:00:00           # Run time (hh:mm:ss)
#SBATCH -p batch              # Desired partition

# Remark: can use -C intel to require the 'intel' feature

# Note that the output from each step will be saved to a unique
# file: %J maps to jobid.stepid

rm -rf build mod bin
mkdir build
cd build
module load intel-oneapi-compilers intel-oneapi-mkl netlib-lapack hpcx-mpi
export FC=ifx
export OMPI_FC=${FC}
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_OMP=ON ..
make -j
cd ../bin
export OMP_PROC_BIND=close
export OMP_PLACES=cores
srun ./main.exe
echo "YUEWU: All Steps completed."
