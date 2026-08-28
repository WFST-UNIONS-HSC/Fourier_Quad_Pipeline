# Fourier_Quad Pipeline

> English: [README.md](README.md)

## 概述

基于Fourier_Quad 剪切测量方法的 Fourier_Quad 弱透镜 MPI 软件。本仓库提供当前 C++17 与旧版
Fortran 两种实现，每种实现均有 Standard 和 Lite 版本。

## Pipeline 版本

| Pipeline | 源码目录 | 用途 |
|---|---|---|
| C++ Lite | [`cpp_Lite`](cpp_Lite/) | 删除未使用替代分支的 C++17 生产路径。 |
| C++ Standard | [`cpp_Standard`](cpp_Standard/) | 保留可选科学路径的完整 C++17 分支集合。 |
| F77 Lite | [`f77_Lite`](f77_Lite/) | 精简的旧版 Fortran 生产路径。 |
| F77 Standard | [`f77`](f77/) | 完整的旧版 Fortran 分支集合。 |

四个版本都生成名为 `Fourier_Quad_Pipe` 的可执行文件。

## 快速开始

### 1. 下载 Release 源码包

打开 [GitHub Releases](https://github.com/WFST-UNIONS-HSC/Fourier_Quad_Pipeline/releases)，
只下载准备运行的 Pipeline 对应源码包：

| Pipeline | Release 源码包 |
|---|---|
| C++ Lite | `cpp_Lite.zip` |
| C++ Standard | `cpp_Standard.zip` |
| F77 Lite | `f77_Lite.zip` |
| F77 Standard | `f77.zip` |

推荐使用固定 Release，以便记录并复现实验实际使用的软件版本。

### 2. 准备运行环境

根据所选版本安装下文 [C++](#c-pipeline) 或 [F77](#f77-pipeline) 环境表中的依赖。

### 3. 准备输入数据

准备[输入数据要求](#输入数据要求)中的四类数据。前三类为必选；DQ masks 是否需要
取决于所选 Pipeline 配置。

### 4. 配置 Pipeline

- C++ 用户：检查所选源码包的 `config/*.hpp`。完整 Standard/Lite 参数对照见
  [CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md)。
- F77 用户：按照 [F77 指南](F77_GUIDE_CN.md)修改 `para.inc`、`cust_para.inc` 和
  `sig_para.inc`。

### 5. 编译和运行

进入解压后的源码包目录，最短入口为：

```bash
# C++ Lite 或 C++ Standard
make -j4
./Fourier_Quad_Pipe --help

# F77 Lite 或 F77 Standard
make LAPACK_LIB_DIR=/path/to/lapack/lib \
     CFITSIO_LIB_DIR=/path/to/cfitsio/lib
mpirun -np 4 ./Fourier_Quad_Pipe /data/work/expo_gband.list
```

完整参数、阶段选择和运行示例见 [C++ 指南](CPP_GUIDE_CN.md)或
[F77 指南](F77_GUIDE_CN.md)。

## 运行环境前置条件

### C++ Pipeline

| 类别 | 前置条件 | 普通 Linux 运行 | HPC / Slurm 运行 | 说明 |
|---|---|---:|---:|---|
| 操作系统 | 64-bit Linux | 必需 | 必需 | 推荐现代 Linux 发行版。 |
| C++ 编译器 | 支持 C++17 的 MPI C++ wrapper | 必需 | 必需 | 当前容器工具链使用 GCC 12.3.0。 |
| MPI | OpenMPI 或兼容 MPI 实现 | 必需 | 必需 | 当前容器使用 OpenMPI 4.1.8。 |
| CFITSIO | CFITSIO 开发库 | 必需 | 必需 | FITS/FZ 读写；容器版本为 4.6.4。 |
| FFTW3 | 双精度与单精度 FFTW3 库 | 必需 | 必需 | Fourier 变换；容器版本为 3.3.11。 |
| Eigen3 | Eigen3 头文件 | 必需 | 必需 | C++ 数值运算；容器版本为 3.4.0。 |
| BLAS / LAPACK | BLAS 与 LAPACK | 必需 | 必需 | 线性代数。 |
| 编译工具 | `make`、shell、标准 GNU 工具 | 必需 | 必需 | 编译与辅助脚本。 |
| 共享文件系统 | 所有 rank/节点可见相同输入输出 | 否 | 必需 | 多节点任务必需。 |
| Slurm | 站点调度器与 `srun` | 否 | Slurm 任务必需 | 批处理分配与启动。 |
| PMI2 兼容启动 | Slurm `pmi2` 与兼容 MPI | 否 | 使用仓库 HPC 容器路径时必需 | 用于 Slurm 直接启动 MPI。 |
| Apptainer / Singularity | 无 root 容器运行时 | 否 | 使用仓库 HPC 容器路径时必需 | 计算节点必须可用。 |
| 可写处理目录 | 所有 rank 可写的共享输出/工作目录 | 必需 | 必需 | 原始数据与输出不应混放。 |

### F77 Pipeline

| 类别 | 前置条件 | 普通 Linux 运行 | HPC / Slurm 运行 | 说明 |
|---|---|---:|---:|---|
| 操作系统 | 64-bit Linux | 必需 | 必需 | F77 Pipeline 的目标平台。 |
| Fortran 编译器 | 兼容的 GNU Fortran 工具链与 `mpif77` | 必需 | 必需 | 固定容器使用 GNU Fortran 4.8.5。 |
| MPI | MPICH 或兼容的 Fortran MPI 实现 | 必需 | 必需 | 固定容器使用 MPICH 4.1.2。 |
| CFITSIO | CFITSIO 库 | 必需 | 必需 | FITS/FZ 读写；容器版本为 4.3.1。 |
| Fourier 例程 | 仓库自带 `FFTPACK.f` | 已包含 | 已包含 | 当前 F77 Makefile 不要求外部 FFTW 库。 |
| BLAS / LAPACK | BLAS 与 LAPACK | 必需 | 必需 | 容器提供 LAPACK 3.8.0 及其 BLAS。 |
| 编译工具 | `make`、shell、标准 GNU 工具 | 必需 | 必需 | 编译与运行脚本。 |
| 共享文件系统 | 所有 task/节点可见相同输入输出 | 否 | 必需 | 多节点任务必需。 |
| Slurm | 站点调度器 | 否 | Slurm 任务必需 | 批处理分配与启动。 |
| MPI 兼容启动 | 与所选容器/宿主 MPI 模式兼容的 launcher | 否 | 必需 | 仓库 runner 支持其文档所述 MPICH 与 PMI 兼容模式。 |
| Apptainer / Singularity | 无 root 容器运行时 | 否 | 使用仓库 HPC 容器路径时必需 | 计算节点必须可用。 |
| 可写处理目录 | 所有 task 可写的共享输出/工作目录 | 必需 | 必需 | 不要把产物写入原始归档。 |

## 输入数据要求

| 输入 | 必选 | 用途 | 最低要求 |
|---|---:|---|---|
| Science images | 是 | Pipeline 对其执行源检测、形状测量和弱透镜处理的科学曝光图像。 | 代码支持的 Science FITS/FZ 数据，且文件组织符合配置的曝光与 CCD 识别规则。 |
| Gaia catalog | 是 | 为天体匹配和测天标定提供高精度天球坐标参考。 | 覆盖 Science images 天区；每个文件第一行为表头，后续行两个字段为数值型 `ra`、`dec`。 |
| External source catalog | 是 | 为外部源匹配和下游处理提供天球位置、`zp` 与光度。 | 包含 `ra`、`dec`、`zp` 和至少一个观测波段的星等。 |
| DQ masks | 取决于配置（输入类别可选） | 标记坏像素、饱和、探测器缺陷及其他无效像素。 | 仅当所选配置不读取 DQ masks 时才可省略。 |

### Science images

Science images 是主要科学曝光图像，不是标定星表或输出目录。它们必须采用所选版本
支持的 FITS/FZ 格式，能够被其 FITS 逻辑读取，并符合详细指南中的曝光/CCD 命名和
列表约定。

### Gaia catalog

Gaia catalog 为天体匹配和测天标定提供精确 RA/Dec 参考位置。它必须覆盖实际 Science
images 天区，并直接存放在配置的 `ASTROMETRY_CAT` 目录下。读取器会跳过第一行表头，
再把后续每行的前两个数值字段读作 RA 和 Dec；数据行可用逗号或空白分隔，额外字段会
被忽略。

**文件名规范：** 
-  `|Dec| < 80°` : `gaia_<p|m><D>_<RR>.cat`
其中`D = floor(|Dec| / 10) + 1`（1-8），`RR = floor(RA / 10)`（00-35，不足两位
补零）。
-  `|Dec| >= 80°` : 不带 RA 后缀 `gaia_<p|m>9.cat`。
- 其中, `p`表示 Dec 非负，`m` 表示 Dec 为负。
- 单个文件必须包含一行表头
> 示例：
> 1. gaia_p1_00.cat 覆盖 `0° <= RA < 10°`、`0° <= Dec < 10°`
> 2. gaia_m3_12.cat 覆盖 `120° <= RA < 130°`、`-30° < Dec <= -20°` 
> 3. gaia_p9.cat 覆盖 `0° <= RA < 360°`、`80° <= Dec <= 90°`
*为了保证位于10°×10°格点边缘的曝光仍能选到星，建议单个星表的dec覆盖上下限在上述基础上
增减2°，ra上下限在0°/30°/60°范围分别增减2°/4°/6°。*

### External source catalog

最低 schema 不绑定固定 survey 或波段：

| 字段 | 含义 |
|---|---|
| `ra` | Right Ascension。 |
| `dec` | Declination。 |
| `zp` | 所选 Pipeline 配置实际使用的 catalog `zp` 量。 |
| 任一观测波段星等 | 所选分析使用的任意一个波段 magnitude。 |

用户可以保留额外颜色、redshift、object class、shape 和 flag，但它们不属于最低输入
要求。实际列位置、列名、delimiter、header 处理和投影由 `config/ExtCatConfig.hpp`
或其已记录的 CLI 覆盖项配置。

**文件名规范：** 
- 1° × 1° 分片 : `des_y6_RA_<RA0>_<RA1>_Dec_<Dec0>_<Dec1>.dat`
- RA 边界固定为三位数字；Dec 边界由 `p` 或 `m` 加两位绝对值组成；每个上边界比下边界大 1 度。
- 单个文件必须包含一行表头
> 示例：
> `des_y6_RA_123_124_Dec_m05_m04.dat` 覆盖 `123° <= RA < 124°`、`-5° <= Dec < -4°`。

### DQ masks（可选）

DQ masks 标记不应参与科学测量的坏像素、饱和、探测器缺陷及其他无效区域。只有当
所选配置关闭 DQ 访问时才可省略。Standard 用户可用 `include_Mask` 选择 mask 模式；
当前 C++ Lite 固定使用逐 CCD DQ masks，因此 Lite 运行仍须提供它们。如果省略 DQ
masks，必须确认所有有效分支和配置路径都不再读取它们。

## C++ 配置参数参考

[CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md) 提供每个配置头文件、
Standard/Lite 默认值、CLI 覆盖方式和重编译要求的完整对照。

## 详细指南

- [C++ Pipeline 指南](CPP_GUIDE_CN.md)
- [F77 Pipeline 指南](F77_GUIDE_CN.md)
- [External source catalog 工具](gen_src_cat/README.md)

## 许可证

仓库自有代码采用 [MIT License](LICENSE)。容器依赖保留其上游许可，详见各容器目录的
第三方声明。
