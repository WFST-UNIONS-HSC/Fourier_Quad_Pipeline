# Fourier_Quad Fortran Slurm/Apptainer Runner

本 runner 将可写 Fortran 源码和数据 bind 到只读 SIF，在批处理节点编译一次，再为每个
MPI rank 启动一个容器命令。

> English: [README.md](README.md)

## 1. 检查站点

```bash
bash inspect-cluster-mpi.sh
```

保存输出，在核对 MPI 构建选项、`srun --mpi=list`、Apptainer 与计算节点网络后再选择
启动方式：

- `mpiexec` 模式要求宿主 MPICH launcher 与 SIF 内 MPICH 4.1.2 应用兼容；
- `srun` 模式要求 Slurm PMI 与 SIF 经过现场验证。pilogin 示例使用
  `srun --mpi=pmi2`；绝不能用宿主 OpenMPI `mpirun` 启动 MPICH 应用。

## 2. 配置共享路径

```bash
cp f77pipeline.env.example f77pipeline.env
```

`f77pipeline.pilogin-openmpi.env.example` 只是站点示例，不是通用默认。设置 SIF、源码、
星表、标定、处理、曝光表、scratch 和可选 cache/runtime 路径；所有宿主路径必须在每个
分配节点同位置可见。

env 文件是 Bash，`HPC_MODULES` 必须保持数组。容器内星表目标要与 `para.inc` 编译路径
一致。

挂载源码的 Makefile 必须在 SIF 内正确解析库。本仓库默认库目录是站点值；提交前应在
私有副本中改为 `/opt/f77stack/lib`，或用其他方式确保镜像内裸
`make -C "$F77_SOURCE_CONTAINER"` 成功。

## 3. 获得与检查 SIF

生产环境优先设置固定 digest 的 `OCI_IMAGE_URI`，然后：

```bash
bash pull-sif.sh
bash run-apptainer.sh --check
```

脚本拒绝覆盖已有 SIF。登录节点不能拉取时，可在获准的 x86_64 Linux 主机生成 SIF 后
传入共享目录。

## 4. 验证 MPI

调整 Slurm 资源模板后执行：

```bash
sbatch mpi-smoke-test.slurm
```

pilogin PMI2 示例应配合对应 env 与 `mpi-smoke-test-pilogin-openmpi.slurm`。新站点必须
依次通过单 rank、单节点多 rank 和多节点测试。

## 5. 运行流水线

```bash
sbatch f77pipeline.slurm
```

无参数时传入 `F77_EXPO_LIST_CONTAINER`；显式位置参数可覆盖：

```bash
sbatch f77pipeline.slurm /data/DataProcess/another_expo.list
```

作业依次检查 bind、可选清理、单次编译、切换处理目录并启动全部 rank。不要让多个作业
同时清理或编译同一源码副本。

资源、模块和 MPI 模式都是站点输入，应成组调整；其他集群的验证不能证明本地
Slurm/PMI 或高速互连兼容。
