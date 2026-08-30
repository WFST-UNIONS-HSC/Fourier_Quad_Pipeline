# Fourier_Quad C++ 流水线指南

本文描述本仓库的 C++ 程序。它使用 CLI 覆盖与编译默认值；与独立的
`Fourier_Quad_Cpp` 仓库不同，本版本没有 INI 配置层。

> English: [CPP_GUIDE.md](CPP_GUIDE.md)

## 版本与目录

[`cpp_Standard`](cpp_Standard/)（C++ Standard）保留平场、掩膜、简化测天、外部/混合
PSF 与 PCA 等可选分支。

[`cpp_Lite`](cpp_Lite/)（C++ Lite）物理删除这些替代路径，
只保留 Gaia 测天、逐 CCD DQ masks、External source catalog 匹配、去混叠、局域多项式
PSF 和无 PCA 路径。由于 C++ Lite 固定使用逐 CCD DQ 分支，其运行必须提供 DQ masks；
C++ Standard 可以选择不读取 DQ 的配置。

每个版本包含 `main.cpp`、`config/`、`include/`、`src/`、`tests/` 与 Makefile。
曝光表、路径、MPI、调度和通用数值工具位于 `include/general/`、`src/general/`；阶段代码
位于六个 `process_*` 子目录。

## 顶层阶段

程序按固定顺序调用启用的阶段：

| 阶段 | CLI | 作用 |
|---|---|---|
| `process_astrocat` | `--run-astrocat` | 将原始两列 Gaia 星表重分块、去重为一度瓦片。 |
| `process_extcat` | `--run-extcat` | 将原始星表重分块为程序标准格式。 |
| `process_init` | `--run-init` | 提取 Science images 与 DQ-mask CCD 并发布曝光表。 |
| `process_main` | `--run-main` | 执行九个数值阶段。 |
| `process_rearr` | `--run-rearr` | 空间分区 `*_all.cat`。 |
| `process_fd` | `--run-fd` | 执行场畸变剪切检验。 |

`process_astrocat`、`process_extcat` 依次各执行一次，其余阶段按数据集顺序执行；首次失败会停止任务。

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

其他安装位置可传入 `CXX`、`STACK_PREFIX` 与 `EIGEN_INCLUDE`。当前 Makefile 只提供
`all` 与 `clean`。应先用 `./Fourier_Quad_Pipe --help` 检查构建，再针对所改配置运行
代表性阶段或数据集。仓库固定容器软件栈使用 GCC 12.3.0、OpenMPI 4.1.8、
CFITSIO 4.6.4、FFTW3 3.3.11 与 Eigen3 3.4.0。


## 配置与 CLI

C++ 程序的编译默认值与固定科学参数位于所选版本的 `config/` 目录。CLI 只覆盖
`ProcessConfig::RuntimeOptions` 中的阶段、I/O、数据集和外部星表选项；修改大多数
`LensingConfig`、`ProcessRearrConfig` 或 `FDConfig` 参数后必须重新编译。
完整选项以 `./Fourier_Quad_Pipe --help` 为准。

CLI 支持 `--name value` 与 `--name=value`。布尔值支持 `true/false`、`1/0`、
`on/off`。首次显式 `--dataset`、`--contains`、`--extcat-contains` 会替换编译列表，
后续重复项追加。一个裸位置参数可作为 `--expo-list` 的兼容写法。

### 集中路径配置

固定输入/输出路径、流程输出与曝光表名称、重排文件名，以及初始化/处理产物的相对目录，
都只在所选版本的 `config/pathconfig.hpp` 中实际定义。科学参数和解析行为仍保留在各自
领域配置头中，原有带命名空间的符号名称不变。

- 必须保留 `AstroCatConfig::ASTROCAT_OUTPUT_DIRECTORY` 从
  `LensingConfig::ASTROMETRY_CAT` 初始化的关系，以及这两个相互耦合但独立的符号。CLI
  `--astrocat-output` 只覆盖运行期生产者目录，不会反向更新 Stage 1 消费的
  `ASTROMETRY_CAT`。
- 必须保留 `ExtCatConfig::EXTCAT_OUTPUT_DIRECTORY` 从
  `LensingConfig::SOURCE_CAT_DEFAULT` 初始化的关系。它的运行期副本同时作为
  `process_extcat` 的输出目录
  和 `process_main` 的输入目录，`--extcat-output` 覆盖该副本。
- `FLAT_PATH` 与 `PSF_PATH` 只存在于 Standard；Lite 已物理删除对应可选分支。
- `ProcessRearrConfig` 的固定文件名和两组 `OutputLayout` 目录数组没有 CLI 覆盖。
  直接修改它们或 `pathconfig.hpp` 中其他默认值后，必须执行 `make clean && make`。

### 常用及随图像数据源变化的参数

下表是运行前应主动检查的参数入口。表中“运行时”表示可用所列 CLI 修改且无需重编译；
“编译时”表示需要修改所选版本的文件并重新执行 `make`。派生尺寸和列号不要单独修改。

| 类别 | 参数（当前默认） | 修改方式 | 何时修改与约束 |
|---|---|---|---|
| 顶层阶段 | `RUN_PROCESS_ASTROCAT`、`RUN_PROCESS_EXTCAT`、`RUN_PROCESS_INIT`、`RUN_PROCESS_MAIN`、`RUN_PROCESS_REARR`、`RUN_PROCESS_FD` | `config/ProcessConfig.hpp`；运行时 `--run-astrocat`、`--run-extcat`、`--run-init`、`--run-main`、`--run-rearr`、`--run-fd` | 选择本次执行的阶段。Standard 默认 `false/false/true/true/true/true`，Lite 默认 `false/false/true/true/false/false`。 |
| Science/DQ 归档与数据集 | `SCIENCE_ROOT`、`DQ_ROOT`、`OUTPUT_ROOT`、`DATASETS`、`CONTAINS` | 路径在 `config/pathconfig.hpp`；数据集/token 在 `config/InitConfig.hpp`；运行时 `--science-root`、`--dq-root`、`--output-root`、`--dataset`、`--contains` | 更换观测归档、文件名前缀、筛选 token、输出根目录时修改。Lite 必须提供逐 CCD DQ masks。 |
| 曝光表与阶段输出 | `EXPO_LIST`、`REARR_OUTPUT_DIRECTORY`、`REARR_OUTPUT_BASE_DIRECTORY`、`REARRANGED_EXPO_LIST_FILENAME`、`REARRANGED_EXPO_LIST_DIRECTORY`、`FD_EXPO_LIST`、`FD_OUTPUT_DIRECTORY`、`FD_OUTPUT_BASE_DIRECTORY` | `config/pathconfig.hpp`；运行时 `--expo-list`、`--rearr-output-dir`、`--rearr-output-base`、`--rearr-list-name`、`--rearr-list-dir`、`--fd-expo-list`、`--fd-output-dir`、`--fd-output-base` | 下游单独运行，或改变重排/FD 输出目录、曝光表位置时修改。 |
| 固定生成布局 | `SKIP_DIRECTORY_NAME`、`SUBCAT_PREFIX`、`SUBCAT_EXTENSION`、`SUMMARY_FILENAME`、`NON_CHIP_BASE_DIRECTORIES`、`CHIP_PRODUCT_DIRECTORIES` | `config/pathconfig.hpp`，编译时 | 仅在发布星表命名或相对输出目录约定变化时修改；重编译并重新生成受影响产物。 |
| Gaia 星表分块 | `ASTROCAT_INPUT_DIRECTORY`、`ASTROCAT_OUTPUT_DIRECTORY`、`ASTROCAT_ADD_HEADER=true`、`ASTROCAT_EXISTING_POLICY=fail` | 路径在 `config/pathconfig.hpp`；行为在 `config/AstroCatConfig.hpp`；运行时 `--astrocat-input`、`--astrocat-output`、`--astrocat-add-header`、`--astrocat-existing` | 更换 Gaia 原始星表或重跑策略时修改。输出选项只控制 `process_astrocat`，不与 `ASTROMETRY_CAT` 校验，也不会传播给它。 |
| Gaia 星表布局 | `AstroCatType=1` | `config/LensingConfig.hpp`，编译时 | `1` 读取旧式大 `gaia_*.cat` 瓦片；`2` 累积读取 `process_astrocat` 生成的一度 `des_y6_*.dat` 瓦片。Stage 1 消费的目录仍应单独写入 `ASTROMETRY_CAT`；修改类型后必须重编译。 |
| 外部星表发现与解析 | `EXTCAT_INPUT_DIRECTORY`、`EXTCAT_OUTPUT_DIRECTORY` | 路径在 `config/pathconfig.hpp`；解析设置在 `config/ExtCatConfig.hpp`；运行时 `--extcat-input`、`--extcat-output` | 更换外部星表文件组织时修改。输出目录不能等于或位于输入目录内。 |
| 外部星表 schema | `EXTCAT_TOTAL_COLUMNS`、`EXTCAT_INPUT_COLUMNS_ONE_BASED`、`EXTCAT_RA_COLUMN_ONE_BASED`、`EXTCAT_DEC_COLUMN_ONE_BASED`、`EXTCAT_ZP_COLUMN_ONE_BASED` | `config/ExtCatConfig.hpp`；投影和 RA/Dec/ZP 列可用 `--extcat-columns`、`--extcat-ra-column`、`--extcat-dec-column`、`--extcat-zp-column` 运行时修改 | 更换 survey 或列顺序时修改。显式投影必须保留 RA、Dec、ZP 和启用阶段消费的字段；改变总列数还需同步审查重排与 FD 列号。 |
| Gaia、外部星表与标定路径 | `ASTROMETRY_CAT`、`SOURCE_CAT_DEFAULT`（有效 `SOURCE_CAT`）、`FLAT_PATH`、`PSF_PATH` | `config/pathconfig.hpp`；`--extcat-output` 可在运行时设置有效 `SOURCE_CAT`，其余为编译时 | 更换 Gaia 瓦片、规范化源星表、平场或外部 PSF 数据源时修改；`--astrocat-output` 与 `ASTROMETRY_CAT` 相互独立；容器内路径必须与 bind 目标一致。 |
| Standard 分支选择 | `ASTROMETRY_trivial=0`、`include_FLAT=0`、`include_Mask=2`、`ext_cat=1`、`ext_PSF=0`、`PSF_type=1`、`PSF_Ms=0` | `config/LensingConfig.hpp`，编译时 | 只有 Standard 可切换这些分支。Lite 已固定为 Gaia、无平场、逐 CCD DQ、外部源星表、帧内 PSF、局域多项式且无 PCA。 |
| 图像与探测器几何 | `npx=3000`、`npy=5000`、`CCD_split=2`、`chipnx=2046`、`chipny=4094`、`pixel_size=0.2628`、`NMAX_CHIP=62`、`NMAX_EXPO=25000` | `config/LensingConfig.hpp`，编译时 | 更换相机、CCD 尺寸、放大器布局、像元尺度或单批曝光规模时修改；几何量必须成组核对。 |
| 数值阶段 | `PROCESS_stage=223092870` | `config/LensingConfig.hpp`，编译时 | 用素因数选择九个主流程阶段；阶段 9（23）必须与阶段 8（19）同时启用。 |
| 源检测与像素阈值 | `saturation_thresh=25000` | `config/LensingConfig.hpp`，编译时 | 更换图像源后，以代表性数据重新标定。 |
| FD 星表列布局 | `col_flags_*`、`col_cra/cdec`、`col_mag_*`、`col_zp`、`col_expo`、`col_ccd` 及派生 `col_*` | `config/FDConfig.hpp`，编译时并需协调 reader/writer | 默认绑定 18 列 DES schema，其后依次为 `EXPO_NUM` 与 `ccD_NUM`。改变外部字段宽度或顺序时，必须同步审查 `ExtCatConfig`、外部星表 reader、重排布局和 FD reader，不能只改一个列号。 |
| FD 探测器规则 | `bad_ccds={2,31,53,61}`、`chip_xmin=50`、`chip_xmax=1990`、`chip_ymin=100`、`chip_ymax=3990` | `config/FDConfig.hpp`，编译时 | 更换相机、坏 CCD 清单或边缘 mask 策略时修改，并同步 `n_bad_ccds`。 |

每个独立参数的 Standard/Lite 默认值、合法值、CLI 覆盖和重编译要求见
[CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md)。修改高耦合参数时，应保留基准
配置，并先用最小代表性数据验证。


## 运行示例

仅运行主流程：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-main true \
  --expo-list /data/work/expo_gband.list
```

初始化并运行主流程：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-init true --run-main true \
  --science-root /data/archive/science --dq-root /data/archive/dqmask \
  --output-root /data/work --dataset g2019:c4d_19 --existing resume
```

仅重分块原始 Gaia 星表：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-astrocat true --run-extcat false --run-init false --run-main false \
  --run-rearr false --run-fd false \
  --astrocat-input /data/raw_gaia --astrocat-output /data/gaia/tiles \
  --astrocat-add-header true --astrocat-existing fail
```

上述 `--astrocat-output` 只控制该阶段；如需 `process_main` 消费此目录，应另行设置
`LensingConfig::ASTROMETRY_CAT`。

仅生成适用于本程序的规范化外部源星表：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe \
  --run-extcat true --run-init false --run-main false \
  --run-rearr false --run-fd false \
  --extcat-input /data/raw_catalogs --extcat-output /data/catalogs/tiles
```

至少启用一个顶层阶段。初始化成功后，后续阶段自动使用生成的绝对路径
`expo_<target>.list`。多个数据集按顺序独立执行。

## 输入与输出

请先按照[顶层输入数据要求](README_CN.md#输入数据要求)准备 Science images、Gaia
catalog、External source catalog 与取决于配置的 DQ masks。顶层章节是唯一的最低 schema
约定；本节只说明 C++ 运行布局和产物。

曝光表每个非空记录包含一个 CCD 列表路径，可带兼容 CCD 数量。初始化器原地读取归档，
在每个数据集下创建 `science/`、`dqmask/`、`stamps/`、`result/`，并发布曝光/fits
列表与 manifest。

主要结果为：

```text
<dataset>/result/<exposure>_all.cat
<dataset>/<rearr-output-dir>/subcat_*.cat
<dataset>/<fd-output-dir>/FD_test_comb.dat
```

`_all.cat`是按曝光为单位的剪切目录，默认包含外部星表字段、原始 1-based
`EXPO_NUM`、1 个 CCD 编号和 25 个流水线字段，共 45 列。schema 升级后需重新生成
Stage 9、rearr 与 FD 产物；旧 44 列数据不能与新版混用。
`subcat_*.cat`是按RA/DEC 重新分块的目录，单个源的所有测量记录连续排列，便于快速去重。
`FD_test_comb.dat`是程序process FD生成的场畸变测试表格文件，用于矫正剪切测量。

常见错误包括：
阶段 9 未同时启用阶段 8；
External source catalog 输出位于输入目录内；
投影遗漏已启用阶段消费的字段；
在 C++ Lite 中启用已删除分支；
在容器参数中使用宿主路径。

参数和输出列见 [CPP_PIPELINE_PARAMETERS.md](CPP_PIPELINE_PARAMETERS.md)。容器见
[cpp_docker/README-CN.md](cpp_docker/README-CN.md)，Slurm 见
[runner 中文指南](cpp_docker/runner/README-CN.md)。
