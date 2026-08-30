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

### 常改 runner 参数

`cpppipeline.env` 中的宿主路径必须在所有分配节点同位置可见；成对的容器路径必须与
程序配置或 CLI 使用的路径一致。

| 参数或位置 | 通常如何修改 | 约束 |
|---|---|---|
| `OCI_IMAGE_URI` / `CPP_DOCKER_ARCHIVE` | 按 OCI 拉取或本地 Docker archive 构建二选一设置。 | `OCI_IMAGE_URI` 指向远程镜像；archive 必须是计算环境可见的 x86_64 镜像归档。 |
| `CPP_SIF`、`CPP_SIF_SHA256_EXPECTED` | 设为共享 SIF 路径和可选预期 SHA256。 | SIF 路径在所有节点一致；生成后再填写校验值。 |
| `CPP_SOURCE_HOST/CONTAINER` | 指向 `cpp_Standard` 或 `cpp_Lite` 及容器内 `/workspace/src_pipe`。 | 宿主源码需可写以保存编译产物。 |
| `SCIENCE_ROOT_*`、`DQ_ROOT_*` | 按初始化阶段使用的 Science/DQ 归档设置。 | 只在需要时 bind；Lite 的主流程需要逐 CCD DQ。 |
| `ASTROMETRY_CAT_*`、`SOURCE_CAT_*`、`FLAT_PATH_*` | 指向 Gaia、External source catalog 和平场目录。 | 容器路径必须与编译路径或受支持的 CLI 覆盖一致。 |
| `PROCESS_DATA_*`、`CPP_EXPO_LIST_CONTAINER` | 设为共享可写处理目录和默认曝光表容器路径。 | 曝光表必须位于已 bind 的容器路径。 |
| `EXTCAT_INPUT_*`、`REARR_OUTPUT_*`、`EXPOLIST_DIR_*`、`FD_OUTPUT_*` | 只为需要独立目录的阶段设置。 | 未设置时不 bind；优先使用 `PROCESS_DATA` 下默认位置。 |
| `HPC_SHARED_SCRATCH_HOST`、`APPTAINER_CACHE_DIR`、`APPTAINER_TMP_DIR` | 设为站点共享可写 scratch/cache/tmp。 | 必须有空间且计算节点可访问。 |
| `APPTAINER_BIN`、`HPC_MODULES`、`SITE_ENV_SCRIPT` | 按站点命令和模块环境设置。 | `HPC_MODULES` 必须保持 Bash 数组；不得注入宿主 MPI 库。 |
| `MPI_LAUNCH_MODE=srun`、`SLURM_MPI_TYPE=pmi2`、`SRUN_ARGS=()` | 保持 runner 的启动方式，必要时在 `SRUN_ARGS` 增加站点参数。 | 本 runner 要求 `srun` + `pmi2`；数组语法不能改成普通字符串。 |
| `HPC_EXTRA_BINDS=()`、`HPC_PASSTHROUGH_ENV`、`HPC_CONTAINER_ENV=()` | 仅在额外文件或环境变量确有需要时设置。 | 保留 Bash 数组格式，并最小化传入 `--cleanenv` 容器的宿主状态。 |
| `CPP_BUILD_JOBS`、`CPP_MAKE_CLEAN`、`CPP_EXECUTABLE` | 调整编译并行度、是否先清理和可执行文件容器路径。 | 不要让多个作业同时清理或编译同一源码副本。 |
| `CPP_IMAGE_ID_EXPECTED`、编译器/MPI 预期版本 | 只在更换镜像软件栈时同步修改。 | 必须与镜像内实际版本成组一致。 |
| `cpppipeline.slurm` 中的 `#SBATCH` | 修改 partition、nodes、ntasks、ntasks-per-node、CPU、内存、time 和日志路径。 | 遵守站点策略，并保持任务数与 MPI rank 规划一致。 |

env 文件是 Bash；保留数组写法，并保持 `MPI_LAUNCH_MODE=srun`、
`SLURM_MPI_TYPE=pmi2`。模块可提供 Slurm/Apptainer，但不得把宿主 MPI 库注入应用。

## 获得与验证 SIF

任选一种来源：

```bash
sbatch build-sif.slurm          # 本地已有 Docker archive
bash pull-sif.sh                # 下载远程仓库 OCI 镜像
```

二者都拒绝覆盖已有 SIF，并生成 SHA256 sidecar。随后建议执行：

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

运行 `process_astrocat` 时，在 `HPC_EXTRA_BINDS` 中把原始 Gaia 目录设为只读 bind，
并把生成分片写到可写的 `PROCESS_DATA_CONTAINER` 下（或另一个明确的可写 bind）。
用 `--astrocat-input` 和 `--astrocat-output` 传入对应容器路径。输出参数与编译期
`LensingConfig::ASTROMETRY_CAT` 相互独立；消费作业须另行把生成目录 bind 到该路径、
设置 `LensingConfig::AstroCatType = 2` 并重新编译。
