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
