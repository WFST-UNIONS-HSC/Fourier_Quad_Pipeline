# Fourier_Quad Fortran 流水线指南

`f77` 与 `f77_Lite` 是旧版 MPI 流水线。程序完全由编译期参数配置，只接受一个位置参数
曝光表路径。

> English: [F77_GUIDE.md](F77_GUIDE.md)

## 选择版本

- `f77`（F77 Standard）：完整分支，包含 `00_psf_module.f` 中的可选 PCA/多尺度 PSF 存储。
- `f77_Lite`（F77 Lite）：固定生产路径，已删除替代测天、平场、掩膜、源、PSF、去混叠、混合与
  PCA 分支；其阶段 1 固定为历史 Type-1 随机/排序背景与噪声估计路径。

两者均生成 `Fourier_Quad_Pipe`。

## 数值阶段

`para.inc` 的 `PROCESS_stage` 是阶段素数的乘积：

| 阶段 | 素数 | Fortran 入口 | 内容 |
|---:|---:|---|---|
| 1 | 2 | `pre_process` | 背景/噪声预处理 |
| 2 | 3 | `proc_astrometry` | Gaia 测天 |
| 3 | 5 | `proc_source` | 源检测与恒星候选体 |
| 4 | 7 | `proc_FourierT_st1` | 恒星候选体功率谱 |
| 5 | 11 | `proc_PSF` | PSF 建模 |
| 6 | 13 | `proc_FourierT_st2` | 星系功率谱 |
| 7 | 17 | `proc_shear` | Fourier_Quad 剪切估计量 |
| 8 | 19 | `proc_info` | 曝光统计 |
| 9 | 23 | `proc_comb` | 星表合并与标定 |

默认乘积 `223092870` 启用全部阶段。运行阶段 9 时应保留阶段 8，以便合并星表获得有效
曝光统计。

## 配置

请先按照[顶层输入数据要求](README_CN.md#输入数据要求)准备 Science images、Gaia
catalog、External source catalog 和取决于配置的 DQ masks。F77 Lite 使用其冻结的 mask
路径；除非已确认所选源码和配置不读取 DQ，否则不要省略 DQ masks。

所有设置修改后都必须重编译：

| 文件 | 作用 |
|---|---|
| `para.inc` | 阶段、星表/标定路径、stamp 几何、分支、阈值、容量和星表索引 |
| `cust_para.inc` | CCD 几何与 Standard PCA 参数 |
| `sig_para.inc` | 仅 Standard 使用的稳健 mode-bar 噪声平面估计器；Lite 中不存在 |

至少检查 `PROCESS_stage`、`ASTROMETRY_CAT`、`SOURCE_CAT` 和启用时的
`FLAT_PATH`。容器内这些字符串必须使用容器路径，并与 bind 目标一致。

`SOURCE_CAT_TILE_PREFIX` 控制外部星表 tile 文件名前缀（默认 `extern_`）。两个版本的
`strl` 均为 512；过长的曝光表、CCD 表、输入路径或生成产物路径会直接报错，而不会静默截断。

Lite 的冻结分支写在 `para.inc` 开头；仅添加参数不能恢复已删除代码。Lite 不再包含
`PreprocsType` selector，也不依赖 F6/mode-bar 阶段 1 实现。

## 编译

需要 `mpif77`、CFITSIO、LAPACK、BLAS；Fourier 例程从仓库自带的 `FFTPACK.f`
编译，因此当前 Makefile 不要求外部 FFTW 库。普通 Linux 与 HPC 条件见顶层 F77
环境表。

```bash
cd f77                         # 或 f77_Lite
make LAPACK_LIB_DIR=/path/to/lapack/lib \
     CFITSIO_LIB_DIR=/path/to/cfitsio/lib
```

若库文件不叫 `libcfitsio.so`，再传入完整的 `CFITSIO_LIB`。Makefile 仅有 `all` 与
`clean`。其默认库目录是站点路径，通用 Linux 与容器构建应显式传入覆盖值。

## 曝光表与运行

每条曝光记录必须包含 CCD 列表路径和 CCD 数量：

```text
/data/work/stamps/123456.list 60
/data/work/stamps/123457.list 59
```

曝光表与各 CCD 表中的空行都会被忽略。CCD 表内的 science 路径须遵循
`<dataset>/science/<exposure>/<chip>.fits`，流水线会从路径末尾移除三级来得到
`<dataset>`。曝光数与 CCD 数分别受 `NMAX_EXPO`、`NMAX_CHIP` 边界检查。

运行时只传一个位置参数：

```bash
mpirun -np 8 ./Fourier_Quad_Pipe /data/work/expo_gband.list
```

Fortran 程序没有 `--help`、`--run-*` 或 `--config`；请在 `para.inc` 选择阶段后重编译。

## 输出

各阶段在 CCD 表指向的数据集树下写中间产物。Standard 与 Lite 现在使用与
`cpp_Standard`、`cpp_Lite` 相同且区分大小写的产物布局：

| 产物族 | 位置 |
|---|---|
| DQ mask | `dqmask/<exposure>/<exposure>_<ccd>.fits` |
| 归一化图与测天 | `stamps/Norm/<exposure>/`、`astrometry/dat_Astro/<exposure>/`、`astrometry/Head/`、`astrometry/dat_Chk/` |
| 源与恒星中间产物 | 按产品拆分的 `stamps/cat_Orig`、`stamps/dat_*`、`stamps/fits_*`；逐 CCD 产物包含 `<exposure>/` 子目录 |
| 曝光级 PSF 产物 | `stamps/dat_StarInfo`、`stamps/fits_StarP`、`stamps/fits_PsfSrc`、`stamps/dat_StarComp`、`stamps/dat_Rescale` |
| Standard PCA 产物 | `stamps/dat_StarXY/<exposure>/`、`stamps/fits_PsfResi/<exposure>/`、`stamps/dat_Pcs/`、`stamps/dat_StarCompV2/` |
| 最终星表 | `result/<exposure>_all.cat` |

F77 程序不会创建目录。运行前请使用 init_program 初始化数据集树，或自行创建完整产物目录。
阶段 8 仍在曝光表同目录写父级汇总 `expo_info.dat`。阶段 9 的每一行都包含从 1 开始的
`EXPO_NUM`、物理 `ccD_NUM`、24 列剪切记录和 `Chi2`；启用外部星表时，其字段位于这些列之前。

请使用可写处理目录，不要原地修改原始归档或星表，并确保所有 MPI rank 看到相同绝对路径。

本地容器见 [f77_docker/README-CN.md](f77_docker/README-CN.md)，无 root Slurm 部署见
[runner 中文指南](f77_docker/runner/README-CN.md)。镜像不包含源码和数据，二者均由
宿主 bind 挂载。
