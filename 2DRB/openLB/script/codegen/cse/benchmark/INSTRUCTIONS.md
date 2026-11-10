# Compiling BabelStream

From this directory
```
cd babelStream
```

## HIP
```
cmake -Bbuild -H. -DMODEL=hip -DCMAKE_CXX_COMPILER=hipcc -DCMAKE_CXX_FLAGS="-O3 -std=c++20 --offload-arch=gfx90a -DFLOAT"
cmake --build build
```
## SYCL
```
cmake -Bbuild -H. -DMODEL=sycl2020-acc -DSYCL_COMPILER=ONEAPI-ICPX -DCXX_EXTRA_FLAGS="-O3 -std=c++20 -DFLOAT -fsycl -fsycl-device-code-split=per_kernel -fsycl-targets=intel_gpu_pvc"
cmake --build build
```

# Compiling the Dynamics

From this directory

## LUMI
```
cp lumi/compile.slurm .
sbatch compile.slurm
```

## AURORA
```
cp aurora/compile.pbs .
qsub compile.slurm
```

# Running the Babel Stream and the Dynamics Benchmark

## LUMI
``` 
sbatch benchmark_all.slurm $MILLION_CELLS
```

## AURORA
```
qsub -v MCELLS=$MILLION_CELLS benchmark_all.pbs
``` 
