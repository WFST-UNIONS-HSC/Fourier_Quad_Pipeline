# Fourier_Quad_Pipeline

基于 **Fourier\_Quad** 方法的弱引力透镜剪切测量流水线，提供 **Fortran 77 (旧版)** 和
**C++ (新版)** 两套独立实现，并附带开箱即用的 **Docker** 构建环境与 **HPC**
（Slurm/Apptainer）运行脚本。

> **English documentation**: See [README.md](README.md)

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

*推荐使用适用于现代环境的Cpp版本。*
每套代码库各有两个变体：

- **完整版**（`f77`、`cpp_Standard`）：包含全部功能，含多尺度 / PCA PSF 重建
  （`PSFRecons` / `proc_psfreconsV2.f`）。
- **精简版**（`f77_Lite`、`cpp_Lite`）：冻结分支后的精简版本。将 8 个编译期开关
  固定为典型的生产环境取值，并物理删除所有未选中分支的代码。详见
  `cpp_Lite/REFACTOR_NOTES.md`。

---


## 快速入门

两种运行方式。选择一个代码库（`f77` / `f77_Lite` 或
`cpp_Standard` / `cpp_Lite`）与下列任一模式；完整源码、参数、Docker 与 HPC runner 说明见
[`F77_GUIDE_CN.md`](F77_GUIDE_CN.md) 或 [`CPP_GUIDE_CN.md`](CPP_GUIDE_CN.md)。

首先需要在Release中下载相关源码压缩包，

### 方式 A - 源码编译运行

手动安装安装工具链，编译后直接运行可执行程序。
```bash

# Fortran 77：gfortran + MPICH + CFITSIO + LAPACK/BLAS
cd f77            # 或 f77_Lite
# 编辑 para.inc（路径、PROCESS_stage 等）
make
mpirun -np 4 ./Fourier_Quad_Pipe expo_list.list

# C++17：g++ + MPI + CFITSIO + FFTW3 + LAPACK/BLAS + Eigen3
cd cpp_Standard   # 或 cpp_Lite
# 编辑 include/process_main/LensingConfig.hpp 与 include/ProcessConfig.hpp
make -j4
mpirun -np 4 ./Fourier_Quad_Pipe --expo-list expo_list.list
```

依赖与 Makefile 覆盖变量：[`F77_GUIDE_CN.md`](F77_GUIDE_CN.md) /
[`CPP_GUIDE_CN.md`](CPP_GUIDE_CN.md)。完整 C++ 参数参考见
[`CPP_PIPELINE_PARAMETERS.md`](CPP_PIPELINE_PARAMETERS.md)。

### 方式 B - 下载 runner 并拉取镜像运行

适用于不想手动配置环境、或F77程序配置运行环境不方便的情况。

从 GHCR 拉取已发布镜像（或用 Docker 环境本地构建），再用随附的 `runner/`
脚本在本地或 HPC 集群运行。

```bash
# 拉取已发布镜像（或：cd f77_docker && docker compose build）
docker pull ghcr.io/wfst-unions-hsc/fourier_quad_pipeline/f77pipeline:sha-2d88686
docker pull ghcr.io/wfst-unions-hsc/fourier_quad_pipeline/cpppipeline:sha-2d88686

# HPC（Slurm + Apptainer）：转为 SIF 并提交
bash f77_docker/runner/pull-sif.sh        # 或 cpp_docker/runner/pull-sif.sh
cp f77_docker/runner/f77pipeline.env.example f77_docker/runner/f77pipeline.env
# 按你的集群编辑 .env 路径
sbatch f77_docker/runner/f77pipeline.slurm
```

Docker 环境与 HPC runner 详细说明：[`F77_GUIDE_CN.md`](F77_GUIDE_CN.md) /
[`CPP_GUIDE_CN.md`](CPP_GUIDE_CN.md)。

---

## 贡献
这个项目是张骏教授一系列Fourier Quad方法的开源实现：

- Zhang, J. (2007). Measuring the cosmic shear in Fourier space: Measuring the cosmic shear in Fourier space. Monthly Notices of the Royal Astronomical Society, 383(1), 113–118. https://doi.org/10.1111/j.1365-2966.2007.12585.x
- Zhang, J., Luo, W., & Foucaud, S. (2015). Accurate shear measurement with faint sources. Journal of Cosmology and Astroparticle Physics, 2015(01), 024–024. https://doi.org/10.1088/1475-7516/2015/01/024
- Zhang, J., Zhang, P., & Luo, W. (2017). APPROACHING THE CRAMÉR–RAO BOUND IN WEAK LENSING WITH PDF SYMMETRIZATION. The Astrophysical Journal, 834(1), 8. https://doi.org/10.3847/1538-4357/834/1/8



## 许可证

仓库原创文件采用 [MIT 许可证](LICENSE)。下载的依赖库及 GCC 兼容性补丁仍受其上游
许可证约束，详见：
- [`f77_docker/THIRD_PARTY_NOTICES.md`](f77_docker/THIRD_PARTY_NOTICES.md)
- [`cpp_docker/THIRD_PARTY_NOTICES.md`](cpp_docker/THIRD_PARTY_NOTICES.md)

