#!/bin/bash
#SBATCH --mem=50G  # 请求 50GB 内存
#SBATCH --job-name=letsgo
#SBATCH --ntasks=28
#SBATCH --time=24:00:00   # 最大运行时间
#SBATCH --partition=normal
#SBATCH --output=/scratch/groups/yyanmin/zongmingnew2/NeuronTransportIGA/nsvms_src/letsgo_$
#SBATCH --error=/scratch/groups/yyanmin/zongmingnew2/NeuronTransportIGA/nsvms_src/letsgo_%$


# 加载所需模块
module load openmpi/4.1.2
module load petsc/3.18.5

# 运行程
mpiexec -np 28 ./nsvms /scratch/groups/yyanmin/zongmingnew4/ 28


