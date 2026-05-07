#!/bin/bash
#DSUB -n test
#DSUB -A root.bingxing2.gpuuser102

#DSUB -q root.default
#DSUB -l wuhanG5500
#DSUB --job_type cosched

#DSUB -N 1
#DSUB -R 'cpu=6;gpu=1;mem=45000'
##DSUB -R 'cpu=48;gpu=8;mem=360000'
#DSUB -e %J.out
#DSUB -o %J.out


#=======================加载环境变量=========================#
export MODULEPATH=/home/bingxing2/apps/modulefiles
export MODULEPATH=/home/bingxing2/gpuuser102/apps/modulefiles:$MODULEPATH

export C_INCLUDE_PATH=/home/bingxing2/gpuuser102/apps/sys/usr/include:$C_INCLUDE_PATH
export CPATH=/home/bingxing2/gpuuser102/apps/sys/usr/include:$CPATH
export LIBRARY_PATH=/home/bingxing2/gpuuser102/apps/sys/usr/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=/home/bingxing2/gpuuser102/apps/sys/usr/lib64:$LD_LIBRARY_PATH
export PATH=/home/bingxing2/gpuuser102/apps/sys/usr/bin:$PATH
export PKG_CONFIG_PATH=/home/bingxing2/gpuuser102/apps/sys/usr/lib64/pkgconfig:$PKG_CONFIG_PATH

#export NCCL_DEBUG=INFO
#export NCCL_IB_DISABLE=0
#export NCCL_IB_HCA=mlx5_0:1
#export NCCL_IB_GID_INDEX=3

module load openmpi/4.1.1

if [ "${CCSCHEDULER_ALLOC_FILE}" != "" ]; then
    echo " "
    ls -la ${CCSCHEDULER_ALLOC_FILE}
    echo ------ cat ${CCSCHEDULER_ALLOC_FILE}
    cat ${CCSCHEDULER_ALLOC_FILE}
fi

#====================获取机器文件 HOSTFILE====================#
export HOSTFILE=./hostfile.$$
rm -rf $HOSTFILE
touch $HOSTFILE
ntask=`cat ${CCSCHEDULER_ALLOC_FILE} | awk -v fff="$HOSTFILE" '{}
{
    split($0, a, " ")
    if (length(a[1]) >0 && length(a[3]) >0) {
        print a[1]" slots=4" >> fff
        total_task+=a[2]
    }
}END{print total_task}'`
echo "hostfile $HOSTFILE generated:"
echo "-----------------------"
cat $HOSTFILE

nvidia-smi topo -m
#nvidia-smi
#lsmod | grep nv

#ucx_info -v

#ucx_info -d | grep cuda
#mpirun -np 2 ./a.out
#ompi_info | grep ucx

mpirun -np 1 -x UCX_MEMTYPE_CACHE=n -x UCX_RNDV_THRESH=8192 -x LD_LIBRARY_PATH  -mca pml ucx --mca btl ^vader,tcp,openib,uct -x UCX_IB_GPU_DIRECT_RDMA=1 -x UCX_NET_DEVICES=mlx5_0:1 --mca plm_rsh_agent /opt/batch/agent/tools/dstart ./a.out  2>&1

#mpirun -np 2 --hostfile $HOSTFILE -map-by dist -mca rmaps_dist_device mlx5_0 -x UCX_MEMTYPE_CACHE=n -x UCX_RNDV_THRESH=8192  -x CUDA_VISIBLE_DEVICES=4,5,6,7  -x LD_LIBRARY_PATH  -mca pml ucx --mca btl ^vader,tcp,openib,uct -x UCX_IB_GPU_DIRECT_RDMA=1 -x UCX_NET_DEVICES=mlx5_0:1 --mca plm_rsh_agent /opt/batch/agent/tools/dstart ./a.out  2>&1

