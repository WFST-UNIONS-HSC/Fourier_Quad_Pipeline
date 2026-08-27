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

容器内执行：

```bash
make -C /workspace/src_pipe -j4
/workspace/src_pipe/Fourier_Quad_Pipe --help
mpirun -np 4 /workspace/src_pipe/Fourier_Quad_Pipe \
  --run-init false --run-main true --run-rearr false --run-fd false \
  --expo-list /data/DataProcess/expo_list.list
```

程序参数必须使用容器路径。核心 bind 为源码、测天/源星表、平场和可写处理目录；
Science/DQ 归档以及 extcat/rearr/曝光表/FD 挂载只在对应阶段需要时启用。

编译科学分支使用的星表/标定目标必须与 `config/LensingConfig.hpp` 一致。
`--extcat-output` 可覆盖单次调用的外部源星表瓦片路径。

## Slurm

将同一镜像转换为一个 SIF，然后按 [runner 中文指南](runner/README-CN.md) 操作。
支持的启动链为：

```text
srun --mpi=pmi2 -> run-apptainer.sh -> apptainer exec --cleanenv -> Fourier_Quad_Pipe
```

依赖来源与许可见 [SOURCES.md](SOURCES.md) 和
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
