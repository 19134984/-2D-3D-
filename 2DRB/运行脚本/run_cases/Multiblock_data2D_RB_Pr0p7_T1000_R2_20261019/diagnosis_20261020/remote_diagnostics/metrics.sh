#!/bin/bash
export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/compilers/lib:${LD_LIBRARY_PATH:-}
echo 'FIELD original6'; '/data2/XLLi/Multiblock/data2D/RB/diagnosis_20261020/snapshot_metrics.exe' '/data2/XLLi/Multiblock/data2D/RB/Ra1e6/results/buoyancyCavity2DOpenaccMultiblockSnapshot-1000.340531252631tff-2217780itc.bin'
echo 'FIELD seed6'; '/data2/XLLi/Multiblock/data2D/RB/diagnosis_20261020/snapshot_metrics.exe' '/data2/XLLi/Multiblock/data2D/RB/diagnosis_20261020/Ra6_seed_export/results/buoyancyCavity2DOpenaccMultiblockSnapshot-200.268374625151tff-444000itc.bin'
