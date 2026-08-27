# Fourier_Quad C++ Slurm/Apptainer Runner

本 runner 在提供 `pmi2` 的 x86_64 Slurm 集群上启动预编译 SIF。应用使用 SIF 内的
编译器、OpenMPI 与科学库。

> English: [README.md](README.md)

## 前提与配置

确认 `srun --mpi=list` 包含 `pmi2`、计算节点可用 Apptainer/Singularity，且所有
路径在各节点同位置可见。然后执行：

```bash
bash inspect-cluster-mpi.sh
cp cpppipeline.env.example cpppipeline.env
```

设置镜像/SIF、源码、处理目录、星表、平场、cache、tmp。Science/DQ/extcat/rearr/
曝光表/FD bind 只在所选阶段需要时设置。

env 文件是 Bash；保留数组写法，并保持 `MPI_LAUNCH_MODE=srun`、
`SLURM_MPI_TYPE=pmi2`。模块可提供 Slurm/Apptainer，但不得把宿主 MPI 库注入应用。

## 获得与验证 SIF

任选一种来源：

```bash
sbatch build-sif.slurm          # 已审核 Docker archive
bash pull-sif.sh                # 固定 digest 的 OCI 镜像
```

二者都拒绝覆盖已有 SIF，并生成 SHA256 sidecar。随后执行：

```bash
bash run-apptainer.sh --check
sbatch compile-pipeline.slurm
sbatch --nodes=1 --ntasks=2 --ntasks-per-node=2 mpi-smoke-test.slurm
sbatch mpi-smoke-test.slurm
```

## 运行

在脚本名后传入 C++ CLI：

```bash
sbatch cpppipeline.slurm \
  --run-init false --run-main true --run-rearr false --run-fd false \
  --expo-list /data/DataProcess/expo_list.list
```

没有参数时，runner 将 `CPP_EXPO_LIST_CONTAINER` 作为兼容曝光表位置参数。源码 bind
可写以保存编译产物；不要并发编译同一副本。按站点修改 Slurm 资源，但不要改变
`srun --mpi=pmi2` 启动边界。
