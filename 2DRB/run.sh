#!/bin/bash

# ==========================================
# 3D LBM OpenMP 高性能运行脚本 (24核优化版)
# ==========================================

# 1. 内存安全防线：解除局部大数组的栈空间限制 (防止 Segmentation Fault)
ulimit -s unlimited

# 2. OpenMP 线程设置 (需与物理核数匹配)
export OMP_NUM_THREADS=24

# 3. 核心绑定防线：防止线程在虚拟核间跳跃，锁定 LBM 访存命中率
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# 4. 日志实时输出防线：强制 Fortran 取消块缓冲，确保 tail -f 能实时看到数据
export GFORTRAN_UNBUFFERED_ALL=1

# 5. 正式启动程序
echo "---------------------------------------------------"
echo "Starting LBM simulation at $(date)"
echo "---------------------------------------------------"
./a.out