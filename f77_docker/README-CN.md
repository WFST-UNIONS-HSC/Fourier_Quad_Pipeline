# Fourier_Quad Fortran 容器

本目录构建供 `f77` 或 `f77_Lite` 使用的可复现 x86_64 工具链。源码、星表、标定、
处理数据和输出位于镜像外，通过宿主 bind 挂载。

> English: [README.md](README.md)

## 运行环境

| 组件 | 版本 |
|---|---|
| Rocky Linux | 8.10 |
| GCC / GFortran | 4.8.5 |
| MPICH | 4.1.2（`ch4:ofi`） |
| CFITSIO | 4.3.1 |
| LAPACK / reference BLAS | 3.8.0 |

依赖来源和校验值见 [SOURCES.md](SOURCES.md) 与 `checksums.sha256`。

## 构建与验证

### 拉取GHCR镜像

```bash
docker pull ghcr.io/wfst-unions-hsc/fourier_quad_pipeline/f77pipeline:latest
```

### 下载源码编译

下载Release中最新版f77_docker.zip。

```bash
docker build --platform linux/amd64 --build-arg BUILD_JOBS=4 \
  -t f77pipeline-dev:gnu4.8.5 .
bash scripts/verify-image.sh f77pipeline-dev:gnu4.8.5
```

## 本地使用

```bash
cp .env.example .env
# 设置 F77_SOURCE_HOST 以及星表、标定和处理目录。
docker compose run --rm FourierQuad-F77
```

### 常改 `.env` 参数

先复制 `.env.example`，再按实际宿主目录修改下表。`*_HOST` 是宿主路径，
`*_CONTAINER` 是程序在容器内看到的绝对路径。

| 参数 | 通常如何修改 | 约束 |
|---|---|---|
| `IMAGE_NAME` | 设为准备运行或本地构建的 F77 镜像 tag。 | 必须与 `docker build -t` 或已拉取镜像一致。 |
| `BASE_IMAGE` | 更换基础镜像 registry 或固定 digest。 | 只在重建工具链镜像时修改。 |
| `BUILD_JOBS` | 设为镜像依赖构建的并行任务数。 | 按本机 CPU 和内存调整。 |
| `HOST_UID`、`HOST_GID` | 设为当前宿主用户 UID/GID。 | 输出需要由宿主用户直接读写时修改。 |
| `F77_SOURCE_HOST` | 指向 `f77` 或 `f77_Lite` 源码目录。 | 以读写方式挂载到 `/workspace/f77`，用于保存编译产物。 |
| `ASTROMETRY_CAT_HOST/CONTAINER` | 指向 Gaia 星表目录。 | 容器路径必须与 `para.inc` 中的测天星表路径一致。 |
| `SOURCE_CAT_HOST/CONTAINER` | 指向 External source catalog 目录。 | 容器路径必须与 `para.inc` 中的源星表路径一致。 |
| `FLAT_PATH_HOST/CONTAINER` | 指向平场标定目录。 | 仅启用平场分支时需要，并须与 `para.inc` 一致。 |
| `PROCESS_DATA_HOST/CONTAINER` | 指向可写处理目录。 | 保存曝光表、中间文件和结果；程序参数使用容器路径。 |

容器内使用镜像科学库编译挂载源码：

```bash
make -C /workspace/f77 clean
make -C /workspace/f77 \
  LAPACK_LIB_DIR=/opt/f77stack/lib \
  CFITSIO_LIB_DIR=/opt/f77stack/lib -j4
mpiexec -n 4 /workspace/f77/Fourier_Quad_Pipe \
  /data/DataProcess/expo_list.list
```

仓库 Makefile 的默认库目录是站点路径；通用/容器构建必须按上例覆盖，或修改私有
Makefile 副本。`para.inc` 编译进去的路径必须使用容器路径，并与
`ASTROMETRY_CAT_CONTAINER`、`SOURCE_CAT_CONTAINER` 和启用时的
`FLAT_PATH_CONTAINER` bind 一致。

修改源码或 include 后只需重编译可执行文件；只有工具链、依赖、Dockerfile 或兼容补丁
改变时才重建镜像。

## HPC

Docker Compose 只用于本地。Slurm 上将已发布镜像转换为 SIF，再按
[runner 中文指南](runner/README-CN.md) 操作。通用 runner 支持兼容宿主 MPICH 的
`mpiexec` hybrid 模式，或经现场验证的 Slurm PMI 模式。不要用宿主 OpenMPI
`mpirun` 启动 SIF 内链接 MPICH 的应用。

第三方许可见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
