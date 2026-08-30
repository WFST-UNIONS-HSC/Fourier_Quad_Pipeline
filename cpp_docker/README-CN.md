# Fourier_Quad C++ 容器

本目录构建供 `cpp_Standard` 或 `cpp_Lite` 使用的 x86_64 Linux 工具链镜像。源码、
星表、观测数据与输出位于镜像外，通过宿主 bind 挂载。

> English: [README.md](README.md)

## 运行环境

Rocky Linux 8.10 镜像包含 G++ 12.3.0、OpenMPI 4.1.8（PMI2）、
CFITSIO 4.6.4、FFTW 3.3.11、Eigen 3.4.0、LAPACK 3.11.0、
OpenBLAS 0.3.33。

通用 HPC 基线要求 x86_64、Slurm `pmi2`、Apptainer/Singularity、共享文件系统和
可路由 TCP。其他架构、只提供 PMIx 的站点、其他调度器或 vendor fabric 需另行验证。

## 构建与验证

### 拉取GHCR镜像

```bash
docker pull ghcr.io/wfst-unions-hsc/fourier_quad_pipeline/cpppipeline:latest
```

### 下载源码编译

下载Release中最新版cpp_docker.zip。

```bash
docker build --platform linux/amd64 --target runtime \
  --build-arg BUILD_JOBS=4 \
  -t cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2 .
bash scripts/verify-image.sh cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2
```

## 本地使用

```bash
cp .env.example .env
# 设置 CPP_SOURCE_HOST 和启用阶段需要的宿主路径。
docker compose run --rm FourierQuad-CPP
```

### 常改 `.env` 参数

先复制 `.env.example`，再按实际宿主目录和所选阶段修改下表。`*_HOST` 是宿主路径，
`*_CONTAINER` 是程序在容器内看到的绝对路径。

| 参数 | 通常如何修改 | 约束 |
|---|---|---|
| `IMAGE_NAME` | 设为准备运行或本地构建的镜像 tag。 | 必须与 `docker build -t` 或已拉取镜像一致。 |
| `BUILD_JOBS` | 设为构建镜像时允许的并行任务数。 | 按本机 CPU 和内存调整。 |
| `HOST_UID`、`HOST_GID` | 设为当前宿主用户 UID/GID。 | 输出需要由宿主用户直接读写时修改。 |
| `CPP_SOURCE_HOST` | 指向 `cpp_Standard` 或 `cpp_Lite` 源码目录。 | 以读写方式挂载到 `/workspace/src_pipe`，用于保存编译产物。 |
| `SCIENCE_ROOT_HOST/CONTAINER` | 指向 Science image 归档及其容器路径。 | 仅 `process_init` 需要；CLI 中使用容器路径。 |
| `DQ_ROOT_HOST/CONTAINER` | 指向 DQ mask 归档及其容器路径。 | 启用 DQ 访问时必须设置；Lite 必须提供。 |
| `ASTROMETRY_CAT_HOST/CONTAINER` | 指向 Gaia 星表目录。 | 容器路径必须等于编译的 `ASTROMETRY_CAT`。 |
| `SOURCE_CAT_HOST/CONTAINER` | 指向规范化 External source catalog 目录。 | 容器路径必须等于有效 `SOURCE_CAT`，或与 `--extcat-output` 一致。 |
| `FLAT_PATH_HOST/CONTAINER` | 指向平场标定目录。 | 仅启用平场分支时需要；容器路径必须等于编译的 `FLAT_PATH`。 |
| `PROCESS_DATA_HOST/CONTAINER` | 指向可写处理目录。 | 保存曝光表、中间文件和结果；容器内默认 `/data/DataProcess`。 |
| `EXTCAT_INPUT_*`、`REARR_OUTPUT_*`、`EXPOLIST_DIR_*`、`FD_OUTPUT_*` | 只为需要独立挂载的相应阶段设置。 | 使用时同时加载 `compose.optional.yaml`；未设置时使用处理目录下默认位置。 |

容器内执行：

```bash
make -C /workspace/src_pipe -j4
/workspace/src_pipe/Fourier_Quad_Pipe --help
mpirun -np 4 /workspace/src_pipe/Fourier_Quad_Pipe \
  --run-extcat true --run-init false --run-main true --run-rearr false --run-fd false \
  --expo-list /data/DataProcess/expo_list.list
```

程序参数必须使用容器路径。核心 bind 为源码、测天/源星表、平场和可写处理目录；
Science/DQ 归档以及 extcat/rearr/曝光表/FD 挂载只在对应阶段需要时启用。

编译科学分支使用的星表/标定目标必须与 `config/pathconfig.hpp` 一致。
`--extcat-output` 可覆盖单次调用的外部源星表瓦片路径。

运行 `process_astrocat` 时，应通过合适的只读 bind 暴露原始 Gaia 目录，并用
`--astrocat-input` 传入其容器路径。测天星表 bind 是只读的，因此须用
`--astrocat-output` 把输出指定到可写位置，通常放在 `PROCESS_DATA_CONTAINER` 下。
该参数只控制转换器的写出目录，不会与编译期 `LensingConfig::ASTROMETRY_CAT` 校验，
也不会更新它。后续运行若要消费 Type 2 分片，须另行把结果 bind 到编译的测天目录、
设置 `LensingConfig::AstroCatType = 2` 并重新编译。

## Slurm

将同一镜像转换为一个 SIF，然后按 [runner 中文指南](runner/README-CN.md) 操作。
支持的启动链为：

```text
srun --mpi=pmi2 -> run-apptainer.sh -> apptainer exec --cleanenv -> Fourier_Quad_Pipe
```

依赖来源与许可见 [SOURCES.md](SOURCES.md) 和
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
