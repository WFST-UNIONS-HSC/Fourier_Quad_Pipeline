# Fourier_Quad Pipeline

面向 DECam 数据的 Fourier_Quad 弱透镜 MPI 软件。本仓库同时保存 C++17 流水线、
Fortran 旧版、容器工具链、Slurm runner 与外部源星表生成工具。

> English: [README.md](README.md)

## 选择程序

| 目录 | 用途 |
|---|---|
| [`cpp_Standard`](cpp_Standard/) | 保留可选科学分支的完整 C++17 流水线。 |
| [`cpp_Lite`](cpp_Lite/) | 删除未使用分支的 C++17 生产路径。 |
| [`f77`](f77/) | 完整 Fortran 流水线。 |
| [`f77_Lite`](f77_Lite/) | 精简 Fortran 生产路径。 |
| [`gen_src_cat`](gen_src_cat/) | 独立 MPI 星表重分块工具与 DES Y6 TAP 下载器。 |

四个流水线版本都生成 `Fourier_Quad_Pipe`，应在各自目录内编译运行。

## C++ 快速开始

需要支持 C++17 的 MPI C++ 编译器，以及 CFITSIO、FFTW3、Eigen3、LAPACK、BLAS。

```bash
cd cpp_Lite                    # 或 cpp_Standard
make -j4
./Fourier_Quad_Pipe --help
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init false --run-main true --run-rearr false --run-fd false \
  --expo-list /data/work/expo_gband.list
```

本仓库 C++ 程序由 CLI 与 `config/*.hpp` 默认值配置，不支持 `--config` INI。编译默认
路径不适用于当前环境时，应显式传入顶层阶段与路径选项。

## Fortran 快速开始

需要 `mpif77`、CFITSIO、LAPACK、BLAS；重编译前编辑三份 include。

```bash
cd f77                         # 或 f77_Lite
make LAPACK_LIB_DIR=/path/to/lapack/lib \
     CFITSIO_LIB_DIR=/path/to/cfitsio/lib
mpirun -np 4 ./Fourier_Quad_Pipe /data/work/expo_gband.list
```

Fortran 可执行文件只接受一个位置参数曝光表，不支持 C++ 的阶段 CLI。

## 处理流程

C++ 驱动可执行五个顶层阶段：

```text
process_extcat -> process_init -> process_main -> process_rearr -> process_fd
```

C++ 的 `process_main` 与 Fortran 主程序均包含从预处理到星表合并的九个素数门控数值阶段。

## 容器与 HPC

- [`cpp_docker`](cpp_docker/) 提供 C++ 工具链镜像和 PMI2 Slurm/Apptainer runner。
- [`f77_docker`](f77_docker/) 提供固定 GNU 4.8.5/MPICH 工具链，以及通用与 pilogin
  Slurm 模板。

镜像只含工具链与依赖；源码、星表、观测数据和输出均由宿主 bind 挂载。

## 文档

| 文档 | 内容 |
|---|---|
| [C++ 指南](CPP_GUIDE_CN.md) / [English](CPP_GUIDE.md) | 编译、CLI、阶段、输入与输出 |
| [C++ 参数](CPP_PIPELINE_PARAMETERS.md) | CLI/默认值、编译期参数和默认 44 列星表 |
| [F77 指南](F77_GUIDE_CN.md) / [English](F77_GUIDE.md) | Fortran 配置、编译、运行与输出 |
| [外部星表](gen_src_cat/README.md) | 独立重分块工具与 TAP 下载器 |

## 许可证

仓库自有代码采用 [MIT License](LICENSE)。容器依赖保留上游许可，详见各容器目录的
第三方声明。
