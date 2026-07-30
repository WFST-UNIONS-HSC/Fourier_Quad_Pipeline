# Fourier_Quad_Pipeline

基于 **Fourier\_Quad** 方法的弱引力透镜剪切测量流水线，提供 **Fortran 77** 和
**C++** 两套独立实现，并附带开箱即用的 **Docker** 构建环境与 **HPC**
（Slurm/Apptainer）运行脚本。

> **English documentation**: See [README.md](README.md)

---

## 目录

- [项目概述](#项目概述)
- [仓库结构](#仓库结构)
- [流水线阶段](#流水线阶段)
- [参数配置](#参数配置)
- [源码编译](#源码编译)
- [Docker 环境](#docker-环境)
- [HPC 部署](#hpc-部署)
- [CI/CD 自动化](#cicd-自动化)
- [许可证](#许可证)

---

## 项目概述

Fourier\_Quad 方法通过在傅里叶空间计算四阶矩（即"quad"项）来测量引力剪切。
本仓库包含一条完整的、MPI 并行的数据处理流水线，覆盖从图像预处理、天体测量校准、
源检测、PSF 建模到剪切估计、最终星表合并的全流程。

维护两套独立代码库：

| 代码库 | 语言 | 编译器 | MPI | 关键依赖库 |
|---|---|---|---|---|
| `f77` / `f77_Lite` | Fortran 77（+ Fortran 90 模块） | GCC 4.8.5 (gfortran) | MPICH 4.1.2 | CFITSIO 4.3.1, LAPACK 3.8.0, FFTPACK |
| `cpp_Standard` / `cpp_Lite` | C++17 | G++ 12.3.0 | OpenMPI 4.1.8 | CFITSIO 4.6.4, FFTW 3.3.11, Eigen 3.4.0, LAPACK 3.11.0 |

每套代码库各有两个变体：

- **完整版**（`f77`、`cpp_Standard`）：包含全部功能，含多尺度 / PCA PSF 重建
  （`PSFRecons` / `proc_psfreconsV2.f`）。
- **精简版**（`f77_Lite`、`cpp_Lite`）：冻结分支后的精简版本。将 8 个编译期开关
  固定为生产环境取值，并物理删除所有未选中分支的代码。详见
  `cpp_Lite/REFACTOR_NOTES.md`。

---

## 仓库结构

```
Fourier_Quad_Pipeline/
├── f77/                  完整 Fortran 77 流水线源码
├── f77_Lite/             精简 Fortran 77 流水线（冻结分支）
├── cpp_Standard/         完整 C++17 流水线源码
├── cpp_Lite/             精简 C++17 流水线（冻结分支）
├── f77_docker/           f77 流水线的 Docker 构建环境
├── cpp_docker/           C++ 流水线的 Docker 构建环境
├── .github/workflows/    CI/CD 工作流（GHCR 镜像 + GitHub Release）
├── README.md             英文文档
├── README_CN.md          中文文档（本文件）
├── LICENSE               MIT 许可证
└── .gitignore
```

### 源码目录

#### `f77/` - 完整 Fortran 77 流水线

| 文件 | 说明 |
|---|---|
| `main.f` | 主程序入口。初始化 MPI，读取曝光列表，分发流水线各阶段。 |
| `para.inc` | 主参数文件。控制图像尺寸、阶段选择（`PROCESS_stage`）、星表路径、PSF 阶数、stamp 尺寸、阈值及星表列索引。 |
| `cust_para.inc` | 自定义参数：CCD 几何尺寸与 PCA PSF 分解参数。 |
| `sig_para.inc` | F6 mode-bar 噪声平面估计器参数。 |
| `pre_process.f` | **阶段 1**：平场、掩膜处理、背景/噪声估计。 |
| `proc_astrometry.f` | **阶段 2**：基于 Gaia 参考星表的天体测量校准。 |
| `proc_source.f` | **阶段 3**：源检测与 stamp 提取。 |
| `proc_FFT_st1.f` | **阶段 4**：星系 stamp 的第一阶段傅里叶变换。 |
| `proc_PSF.f` | **阶段 5**：基于恒星 stamp 的 PSF 建模（局域多项式拟合）。 |
| `proc_psfreconsV2.f` | PSF PCA 重建（多尺度，仅 `PSF_Ms=1` 时启用）。 |
| `00_psf_module.f` | Fortran 90 模块，用于全局 PSF PCA 存储。 |
| `proc_FFT_st2.f` | **阶段 6**：第二阶段傅里叶变换。 |
| `proc_shear.f` | **阶段 7**：Fourier\_Quad 剪切估计。 |
| `proc_info.f` | **阶段 8**：单次曝光统计信息收集。 |
| `proc_combine_shear_catalog.f` | **阶段 9**：剪切星表合并与标定。 |
| `astrometry_calib.f` | 天体测量校准工具例程。 |
| `ex_star.f` | 恒星提取与分类。 |
| `rw_fits_image.f` | FITS 图像读写例程（基于 CFITSIO）。 |
| `universal.f` | 通用工具函数。 |
| `univ_imag_proc.f` | 图像处理工具。 |
| `press.f` | Numerical Recipes 例程（排序、统计、插值）。 |
| `mpi_routines.f` | MPI 分发辅助函数（`mpi_distribute`）。 |
| `FFTPACK.f` | FFTPACK 5.1 快速傅里叶变换库。 |
| `Makefile` | 构建文件。使用 `mpif77`，链接 LAPACK、BLAS 和 CFITSIO。 |

#### `f77_Lite/` - 精简 Fortran 77 流水线

文件集与 `f77/` 相同，但不含 `00_psf_module.f`（PCA PSF 重建已移除）。8 个编译期
开关均冻结为生产取值，死代码分支已物理删除。

#### `cpp_Standard/` - 完整 C++17 流水线

| 文件 | 说明 |
|---|---|
| `main.cpp` | 主程序入口，与 `f77/main.f` 一一对应。 |
| `LensingConfig.hpp` | 配置常量（等价于 `para.inc` + `cust_para.inc` + `sig_para.inc`）。 |
| `PreProcess.cpp/.hpp` | **阶段 1**：预处理。 |
| `Astrometry.cpp/.hpp` | **阶段 2**：天体测量校准。 |
| `SourceExtractor.cpp/.hpp` | **阶段 3**：源检测与提取。 |
| `FourierTransformSt1.cpp/.hpp` | **阶段 4**：第一阶段傅里叶变换。 |
| `PSFModel.cpp/.hpp` | **阶段 5**：PSF 建模。 |
| `PSFRecons.cpp/.hpp` | PSF PCA 重建（仅 `PSF_Ms=1` 时启用）。 |
| `FourierTransformSt2.cpp/.hpp` | **阶段 6**：第二阶段傅里叶变换。 |
| `ShearMeasurement.cpp/.hpp` | **阶段 7**：Fourier\_Quad 剪切估计。 |
| `ExposureInfo.cpp/.hpp` | **阶段 8**：单次曝光统计。 |
| `CatalogCombiner.cpp/.hpp` | **阶段 9**：星表合并与标定。 |
| `FitsIO.cpp/.hpp` | FITS 读写例程。 |
| `LinearSolve.cpp/.hpp` | 线性代数工具。 |
| `UniversalUtils.cpp/.hpp` | 通用工具函数。 |
| `ImageProcessing.cpp/.hpp` | 图像处理工具。 |
| `NumericalRecipes.cpp/.hpp` | Numerical Recipes 移植（随机数、排序、插值）。 |
| `MPIScheduler.cpp/.hpp` | MPI 初始化与任务分发。 |
| `ExStar.cpp/.hpp` | 恒星提取与分类。 |
| `Makefile` | 构建文件。使用 `mpicxx`，C++17，链接 CFITSIO、FFTW、LAPACK。 |

#### `cpp_Lite/` - 精简 C++17 流水线

文件集与 `cpp_Standard/` 相同，但不含 `PSFRecons.cpp/.hpp`。详见
`cpp_Lite/REFACTOR_NOTES.md` 了解完整变更记录。

### Docker 目录

#### `f77_docker/`

构建与 pilogin 集群工具链匹配的开发镜像：

| 组件 | 版本 |
|---|---|
| 基础镜像 | Rocky Linux 8.10 (glibc 2.28) |
| GCC / G++ / GFortran | 4.8.5（源码编译） |
| MPICH | 4.1.2 (ch4:ofi 设备，Slurm 启动器) |
| CFITSIO | 4.3.1 |
| LAPACK + 参考 BLAS | 3.8.0 |

关键文件：`Dockerfile`、`compose.yaml`、`.env.example`、
`scripts/verify-image.sh`、`scripts/check-public-repo.sh`、`runner/`（HPC 部署
脚本）、`config/`、`patches/`、`checksums.sha256`、`SOURCES.md`、
`THIRD_PARTY_NOTICES.md`。

#### `cpp_docker/`

构建可移植 HPC 运行时镜像：

| 组件 | 版本 |
|---|---|
| 基础镜像 | Rocky Linux 8.10 |
| G++ | 12.3.0 (conda-forge) |
| OpenMPI | 4.1.8 (Slurm PMI2 直启动) |
| CFITSIO | 4.6.4 |
| FFTW | 3.3.11 |
| Eigen | 3.4.0 |
| LAPACK / OpenBLAS | 3.11.0 / 0.3.33 |

关键文件：`Dockerfile`、`compose.yaml`、`pixi.toml`、`.env.example`、
`scripts/verify-image.sh`、`scripts/check-public-repo.sh`、`runner/`（HPC 部署
脚本）、`SOURCES.md`、`THIRD_PARTY_NOTICES.md`。

---

## 流水线阶段

流水线共 9 个阶段。阶段执行由 `PROCESS_stage` 参数（定义于 `para.inc` /
`LensingConfig.hpp`）控制，该参数为素因数之积。当 `PROCESS_stage` 能被某阶段的
素数整除时，该阶段执行。默认值 `2 * 3 * 5 * 7 * 11 * 13 * 17 * 19 * 23` 启用全部
阶段。

| 阶段 | 素数 | 函数 | 说明 |
|---|---|---|---|
| 1 | 2 | `pre_process` / `PreProcess` | 读取 FITS 图像，平场与掩膜校正，背景噪声估计（F6 mode-bar 估计器）。 |
| 2 | 3 | `proc_astrometry` / `Astrometry` | 基于 Gaia 参考星表的天体测量校准，WCS 拟合。 |
| 3 | 5 | `proc_source` / `SourceExtractor` | 源检测、去混叠、stamp 提取。 |
| 4 | 7 | `proc_FFT_st1` / `FourierTransformSt1` | 星系 stamp 第一阶段傅里叶变换。 |
| 5 | 11 | `proc_PSF` / `PSFModel` | 基于恒星 stamp 的 PSF 建模（局域多项式拟合）。可选 PCA 重建（`PSF_Ms=1`）。 |
| 6 | 13 | `proc_FFT_st2` / `FourierTransformSt2` | 第二阶段傅里叶变换。 |
| 7 | 17 | `proc_shear` / `ShearMeasurement` | 基于四阶傅里叶矩的 Fourier\_Quad 剪切估计。 |
| 8 | 19 | `proc_info` / `ExposureInfo` | 收集单次曝光统计（PSF FWHM、恒星数等）。 |
| 9 | 23 | `proc_combine_shear_catalog` / `CatalogCombiner` | 跨曝光合并剪切星表并应用标定校正。 |

如需跳过某阶段，将 `PROCESS_stage` 除以该阶段的素因数即可。例如设置
`PROCESS_stage = 2 * 3 * 5 * 7 * 11 * 13 * 17 * 19`（省去 23）可跳过星表合并。

---

## 参数配置

### Fortran 77（`para.inc`、`cust_para.inc`、`sig_para.inc`）

所有编译期参数分布于三个 include 文件：

- **`para.inc`** - 图像尺寸（`npx`、`npy`）、阶段控制（`PROCESS_stage`）、星表
  路径、PSF 参数、stamp 尺寸、检测阈值及星表列索引。
- **`cust_para.inc`** - CCD 几何（`chipnx`、`chipny`、`Camera_ccd_num`）及 PCA
  PSF 分解参数。
- **`sig_para.inc`** - F6 mode-bar 噪声平面估计器参数。

> `ASTROMETRY_CAT`、`SOURCE_CAT`、`FLAT_PATH` 等路径必须与 Docker `.env` 文件中
> 定义的容器挂载路径一致。

### C++（`LensingConfig.hpp`）

三个 Fortran include 文件中的全部参数整合到 `LensingConfig.hpp` 的单一
`LensingConfig` 命名空间中。星表列索引已调整为 C++ 的 0 基。

---

## 源码编译

### Fortran 77

前置条件：`gfortran`、`mpif77`（MPICH）、CFITSIO、LAPACK/BLAS。

```bash
cd f77          # 或 f77_Lite
# 编辑 para.inc 设置星表路径和参数
make
# 可执行文件：./Fourier_Quad_Pipe
```

Makefile 支持变量覆盖：

```bash
make LAPACK_LIB_DIR=/path/to/lapack/lib \
     CFITSIO_LIB_DIR=/path/to/cfitsio/lib
```

### C++

前置条件：`g++`（C++17）、`mpicxx`（OpenMPI 或 MPICH）、CFITSIO、FFTW3、
LAPACK/BLAS、Eigen3。

```bash
cd cpp_Standard   # 或 cpp_Lite
# 如需修改星表路径，编辑 LensingConfig.hpp
make -j4
# 可执行文件：./Fourier_Quad_Main
```

`cpp_Standard` 的 Makefile 支持可选的 `STACK_PREFIX`：

```bash
make STACK_PREFIX=/opt/cppstack -j4
```

### 运行流水线

```bash
mpirun -np <N> ./Fourier_Quad_Pipe <EXPO_LIST>    # Fortran
mpirun -np <N> ./Fourier_Quad_Main <EXPO_LIST>    # C++
```

`EXPO_LIST` 是一个文本文件，列出曝光名称和芯片数。每个 MPI 进程被分配一部分
曝光进行处理。

---

## Docker 环境

Docker 环境提供可复现的构建工具链，无需手动安装编译器和科学库。

### 快速开始（f77）

```bash
cd f77_docker
cp .env.example .env
# 编辑 .env：设置 F77_SOURCE_HOST、星表路径和 PROCESS_DATA_HOST

docker compose build
docker compose run --rm FourierQuad-F77
```

进入容器后：

```bash
cd /workspace/f77
make
mpirun -np 4 ./Fourier_Quad_Pipe /data/DataProcess/expo_list.list
```

### 快速开始（C++）

```bash
cd cpp_docker
cp .env.example .env
# 编辑 .env：设置 CPP_SOURCE_HOST、星表路径和 PROCESS_DATA_HOST

docker compose build
docker compose run --rm FourierQuad-CPP
```

进入容器后：

```bash
cd /workspace/cpp_Standard
make -j4
mpirun -np 4 ./Fourier_Quad_Main /data/DataProcess/expo_list.list
```

### 验证镜像

每个 Docker 目录附带验证脚本：

```bash
bash f77_docker/scripts/verify-image.sh f77pipeline-dev:gnu4.8.5
bash cpp_docker/scripts/verify-image.sh cpppipeline-dev:gxx12.3-openmpi4.1.8-pmi2
```

Docker 环境详细文档请参见：
- [`f77_docker/README.md`](f77_docker/README.md) / [`f77_docker/README-CN.md`](f77_docker/README-CN.md)
- [`cpp_docker/README.md`](cpp_docker/README.md) / [`cpp_docker/README-CN.md`](cpp_docker/README-CN.md)

---

## HPC 部署

两套 Docker 环境均包含 `runner/` 目录，提供 Slurm/Apptainer 部署脚本。典型流程：

1. 拉取 GHCR 镜像并转换为 SIF：
   ```bash
   bash f77_docker/runner/pull-sif.sh    # 或 cpp_docker/runner/pull-sif.sh
   ```

2. 配置环境：
   ```bash
   cp f77_docker/runner/f77pipeline.env.example f77_docker/runner/f77pipeline.env
   # 编辑集群路径
   ```

3. 审计 MPI 兼容性（只读）：
   ```bash
   bash f77_docker/runner/inspect-cluster-mpi.sh
   ```

4. 提交流水线：
   ```bash
   sbatch f77_docker/runner/f77pipeline.slurm
   ```

HPC runner 详细文档请参见：
- [`f77_docker/runner/README.md`](f77_docker/runner/README.md) / [`f77_docker/runner/README-CN.md`](f77_docker/runner/README-CN.md)
- [`cpp_docker/runner/README.md`](cpp_docker/runner/README.md) / [`cpp_docker/runner/README-CN.md`](cpp_docker/runner/README-CN.md)

---

## CI/CD 自动化

两个 GitHub Actions 工作流分别自动化镜像发布与 Release 发布。

### Docker 镜像 -> GHCR

**工作流：** [`.github/workflows/docker-images.yml`](.github/workflows/docker-images.yml)

**触发：** 推送 `v*` 标签（例如 `git tag v1.0.0 && git push origin v1.0.0`）
或手动触发。

构建并发布两套镜像到 GitHub Container Registry：

| 镜像 | GHCR 路径 | 标签 |
|---|---|---|
| f77 | `ghcr.io/<owner>/<repo>/f77pipeline` | `gnu4.8.5`、`<版本号>`、`latest`、`<sha>` |
| cpp | `ghcr.io/<owner>/<repo>/cpppipeline` | `gxx12.3-openmpi4.1.8-pmi2`、`<版本号>`、`latest`、`<sha>` |

拉取已发布镜像：

```bash
docker pull ghcr.io/syoong-s/fourier_quad_pipeline/f77pipeline:latest
docker pull ghcr.io/syoong-s/fourier_quad_pipeline/cpppipeline:latest
```

### GitHub Release

**工作流：** [`.github/workflows/release.yml`](.github/workflows/release.yml)

**触发：** 推送 `v*` 标签。

创建 GitHub Release（含自动生成的发布说明），并将仓库下每个顶层内容
目录各自打包为一个 zip 文件作为附件：

```
artifacts/
├── f77.zip
├── f77_Lite.zip
├── f77_docker.zip
├── cpp_Standard.zip
├── cpp_Lite.zip
└── cpp_docker.zip
```

### 发布流程

```bash
# 1. 提交所有更改
git add -A
git commit -m "Prepare release v1.0.0"

# 2. 打标签并推送
git tag v1.0.0
git push origin main --tags
```

推送标签后将同时触发两个工作流。

---

## 许可证

仓库原创文件采用 [MIT 许可证](LICENSE)。下载的依赖库及 GCC 兼容性补丁仍受其上游
许可证约束，详见：
- [`f77_docker/THIRD_PARTY_NOTICES.md`](f77_docker/THIRD_PARTY_NOTICES.md)
- [`cpp_docker/THIRD_PARTY_NOTICES.md`](cpp_docker/THIRD_PARTY_NOTICES.md)
