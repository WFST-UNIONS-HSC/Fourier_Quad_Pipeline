# Fourier_Quad C++ 流水线指南

本文描述本仓库的 C++ 程序。它使用 CLI 覆盖与编译默认值；与独立的
`Fourier_Quad_Cpp` 仓库不同，本版本没有 INI 配置层。

> English: [CPP_GUIDE.md](CPP_GUIDE.md)

## 版本与目录

[`cpp_Standard`](cpp_Standard/) 保留平场、掩膜、简化测天、外部/混合 PSF 与 PCA
等可选分支。[`cpp_Lite`](cpp_Lite/) 物理删除这些替代路径，只保留 Gaia 测天、逐 CCD
DQ 掩膜、外部源、去混叠、局域多项式 PSF 和无 PCA 路径。

每个版本包含 `main.cpp`、`config/`、`include/`、`src/`、`tests/` 与 Makefile。
曝光表、路径、MPI、调度和通用数值工具位于 `include/general/`、`src/general/`；阶段代码
位于五个 `process_*` 子目录。

## 顶层阶段

程序按固定顺序调用启用的阶段：

| 阶段 | CLI | 作用 |
|---|---|---|
| `process_extcat` | `--run-extcat` | 将原始文本星表重分块。 |
| `process_init` | `--run-init` | 提取 Science/DQ CCD 并发布曝光表。 |
| `process_main` | `--run-main` | 执行九个数值阶段。 |
| `process_rearr` | `--run-rearr` | 空间分区 `*_all.cat`。 |
| `process_fd` | `--run-fd` | 执行场畸变剪切检验。 |

`process_extcat` 只执行一次，其余阶段按数据集顺序执行；首次失败会停止任务。

`process_main` 用 `PROCESS_stage` 的素因数选择阶段：

| 阶段 | 素数 | 内容 |
|---:|---:|---|
| 1 | 2 | 预处理与 Gaia 匹配 |
| 2 | 3 | 测天 |
| 3 | 5 | 源检测与恒星候选体 |
| 4 | 7 | 恒星候选体功率谱 |
| 5 | 11 | PSF 建模 |
| 6 | 13 | 星系功率谱 |
| 7 | 17 | Fourier_Quad 剪切估计量 |
| 8 | 19 | 曝光统计 |
| 9 | 23 | 星表合并与标定 |

默认 `223092870` 启用全部阶段；阶段 9 必须与阶段 8 同时启用。

## 编译

需要支持 C++17 的 MPI C++ 编译器，以及 CFITSIO、FFTW3、Eigen3、LAPACK、BLAS。

```bash
cd cpp_Lite                    # 或 cpp_Standard
make -j4
./Fourier_Quad_Pipe --help
```

其他安装位置可传入 `CXX`、`STACK_PREFIX` 与 `EIGEN_INCLUDE`。当前 Makefile 目标包括
`all`、`clean`、`test-general-infrastructure`、`test-psf-star-selection`、
`test-psf-model-state` 与组合目标 `test-stage5`。后者执行两个 Stage-5 专项测试。
本地验证环境为 GCC 15.2.0 的 MPI C++ wrapper，并可用 CFITSIO 4.6.3 与
FFTW3 3.3.10；该环境未安装 Eigen3，因此生产源码编译必须使用包含上述全部依赖的
完整站点软件栈。

Stage 5 使用所有同 CCD FWHM-locus 无序配对计算每个候选体的真实最近形态距离，
曝光阈值则单独由“受限大尺寸 reference 与其他 locus 星”的唯一配对估计。合法的首次
PSF 拟合始终作为回退：raw analytic PRESS 只作为诊断，leverage-standardized PRESS
仅驱动有删除上限的可选删除/重拟合事务。关闭 rejection 不会关闭拟合、leverage、
LOO 或 PRESS 诊断。

## 配置与 CLI

CLI 只覆盖 `ProcessConfig::RuntimeOptions` 表示的字段。`config/LensingConfig.hpp`
中的多数科学参数、rearr 常量和 FD 统计参数仍为编译期设置，修改后必须重编译。

CLI 支持 `--name value` 与 `--name=value`。布尔值支持 `true/false`、`1/0`、
`on/off`。首次显式 `--dataset`、`--contains`、`--extcat-contains` 会替换编译列表，
后续重复项追加。一个裸位置参数可作为 `--expo-list` 的兼容写法。

完整选项以 `./Fourier_Quad_Pipe --help` 为准，主要包括五个 `--run-*`、
`--extcat-*`、初始化器路径/数据集/策略，以及 `--expo-list`、`--rearr-*`、`--fd-*`。
`--extcat-output` 同时改变本次调用的外部源星表路径；其他 lensing 路径和分支仍需编译。

## 运行示例

仅运行主流程：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-init false --run-main true --run-rearr false --run-fd false \
  --expo-list /data/work/expo_gband.list
```

初始化并运行主流程：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-init true --run-main true --run-rearr false --run-fd false \
  --science-root /data/archive/science --dq-root /data/archive/dqmask \
  --output-root /data/work --dataset g2019:c4d_19 --existing resume
```

仅生成外部星表瓦片：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-extcat true --run-init false --run-main false \
  --run-rearr false --run-fd false \
  --extcat-input /data/raw_catalogs --extcat-output /data/catalogs/tiles
```

至少启用一个顶层阶段。初始化成功后，后续阶段自动使用生成的绝对路径
`expo_<target>.list`。多个数据集按顺序执行。

## 输入与输出

曝光表每个非空记录包含一个 CCD 列表路径，可带兼容 CCD 数量。初始化器原地读取归档，
在每个数据集下创建 `science/`、`dqmask/`、`stamps/`、`result/`，并发布曝光/fits
列表与 manifest。

主要结果为：

```text
<dataset>/result/<exposure>_all.cat
<dataset>/<rearr-output-dir>/subcat_*.cat
<dataset>/<rearr-output-dir>/catalog_summary.txt
<dataset>/<fd-output-dir>/FD_test_comb.dat
```

阶段 7 输出截至 WCS parity 的 24 个字段；阶段 9 追加曝光 `chi2`。默认最终行宽为
18 个外部字段 + 1 个 CCD 编号 + 25 个流水线字段 = 44。显式外部列投影只改变外部
前缀宽度，RA、Dec、photo-z 必须仍可用。

常见错误包括：阶段 9 未同时启用阶段 8；外部星表输出位于输入目录内；投影缺少
RA/Dec/photo-z；在 Lite 中启用已删除分支；在容器参数中使用宿主路径。

参数和输出列见 [CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md)。容器见
[cpp_docker/README-CN.md](cpp_docker/README-CN.md)，Slurm 见
[runner 中文指南](cpp_docker/runner/README-CN.md)。
