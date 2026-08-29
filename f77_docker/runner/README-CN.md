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

`f77pipeline.pilogin-openmpi.env.example` 只是站点示例，不是通用默认。

### 常改 runner 参数

所有宿主路径必须在每个分配节点同位置可见；成对的容器路径必须与 `para.inc` 中编译的
路径或程序位置参数一致。

| 参数或位置 | 通常如何修改 | 约束 |
|---|---|---|
| `OCI_IMAGE_URI`、`F77_SIF` | 设为远程 F77 镜像和共享 SIF 目标/现有路径。 | SIF 必须在所有节点同位置可见。 |
| `F77_SOURCE_HOST/CONTAINER` | 指向 `f77` 或 `f77_Lite` 及容器内 `/workspace/f77`。 | 宿主源码需可写以保存编译产物。 |
| `ASTROMETRY_CAT_*`、`SOURCE_CAT_*`、`FLAT_PATH_*` | 指向 Gaia、External source catalog 和平场目录。 | 容器路径必须与 `para.inc` 中编译的路径一致。 |
| `PROCESS_DATA_*`、`F77_EXPO_LIST_CONTAINER` | 设为共享可写处理目录和默认曝光表容器路径。 | 曝光表必须位于已 bind 的容器路径；显式位置参数可覆盖。 |
| `HPC_SHARED_SCRATCH_HOST` | 设为 MPI smoke test 可写的共享目录。 | 所有节点可见，并应有足够空间。 |
| `APPTAINER_BIN`、`HPC_MODULES`、`SITE_ENV_SCRIPT` | 按站点的 Apptainer/Singularity 命令和模块环境设置。 | `HPC_MODULES` 必须保持 Bash 数组。 |
| `MPI_LAUNCH_MODE`、`MPI_LAUNCHER`、`SLURM_MPI_TYPE` | 通用模板默认 `mpiexec`；只有现场验证后才切换 `srun`/`pmi2`。 | launcher 必须与 SIF 内 MPICH 4.1.2 兼容；不得使用宿主 OpenMPI `mpirun`。 |
| `HPC_EXTRA_BINDS` | 需要额外宿主文件/目录时填写逗号分隔 bind 列表。 | 保持模板要求的字符串格式，不要改成 Bash 数组。 |
| `FI_PROVIDER`、`FI_PROVIDER_PATH` | 仅在站点诊断要求指定 libfabric provider 时设置。 | 空值保留镜像 MPICH/libfabric 默认行为。 |
| `HPC_SCRUB_OPENMPI_ENV` | 通常保持 `1`。 | 防止宿主 OpenMPI 传输变量污染 MPICH 容器。 |
| `F77_BUILD_JOBS`、`F77_MAKE_CLEAN`、`F77_EXECUTABLE` | 调整编译并行度、是否先清理和可执行文件容器路径。 | 不要让多个作业同时清理或编译同一源码副本。 |
| `f77pipeline.slurm` 中的 `#SBATCH` | 修改 partition、nodes、ntasks、ntasks-per-node、CPU、内存、time 和日志路径。 | 遵守站点策略，并保持任务数与所选 MPI 模式一致。 |

env 文件是 Bash，`HPC_MODULES` 必须保持数组。容器内星表目标要与 `para.inc` 编译路径
一致。

挂载源码的 Makefile 必须在 SIF 内正确解析库。本仓库默认库目录是站点值；提交前应在
私有副本中改为 `/opt/f77stack/lib`，或用其他方式确保镜像内裸
`make -C "$F77_SOURCE_CONTAINER"` 成功。

## 3. 获得与检查 SIF

设置远程仓库的 `OCI_IMAGE_URI`，然后下载 OCI 镜像：

```bash
bash pull-sif.sh
```

脚本拒绝覆盖已有 SIF。随后建议执行：

```bash
bash run-apptainer.sh --check
```

登录节点不能拉取时，可在获准的 x86_64 Linux 主机生成 SIF 后
传入共享目录。

## 4. 验证 MPI

调整 Slurm 资源模板后执行：

```bash
sbatch mpi-smoke-test.slurm
```

pilogin PMI2 示例应配合对应 env 与 `mpi-smoke-test-pilogin-openmpi.slurm`。建议新站点
依次进行单 rank、单节点多 rank 和多节点测试。

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
