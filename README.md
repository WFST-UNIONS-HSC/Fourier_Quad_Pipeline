# Fourier_Quad Pipeline

MPI software for Fourier_Quad weak-lensing processing of DECam data. This
repository keeps the current C++17 pipeline, the legacy Fortran pipeline,
container toolchains, Slurm runners, and external-catalog generation tools in
one place.

> 中文版：[README_CN.md](README_CN.md)

## Choose a program

| Directory | Purpose |
|---|---|
| [`cpp_Standard`](cpp_Standard/) | Full C++17 pipeline with optional scientific branches. |
| [`cpp_Lite`](cpp_Lite/) | C++17 production path with unused branches removed. |
| [`f77`](f77/) | Full Fortran pipeline. |
| [`f77_Lite`](f77_Lite/) | Reduced Fortran production path. |
| [`gen_src_cat`](gen_src_cat/) | Standalone MPI catalog repartitioner and DES Y6 TAP downloader. |

All four pipeline variants build an executable named `Fourier_Quad_Pipe`.
Build and run them from their own directories.

## C++ quick start

Requires an MPI C++17 compiler, CFITSIO, FFTW3, Eigen3, LAPACK, and BLAS.

```bash
cd cpp_Lite                    # or cpp_Standard
make -j4
./Fourier_Quad_Pipe --help
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init false --run-main true --run-rearr false --run-fd false \
  --expo-list /data/work/expo_gband.list
```

This repository's C++ program is configured by CLI plus defaults in
`config/*.hpp`; it does not accept `--config` INI files. Use explicit phase
switches when the compiled site defaults are not appropriate.

## Fortran quick start

Requires `mpif77`, CFITSIO, LAPACK, and BLAS. Edit the three include files
before rebuilding.

```bash
cd f77                         # or f77_Lite
make LAPACK_LIB_DIR=/path/to/lapack/lib \
     CFITSIO_LIB_DIR=/path/to/cfitsio/lib
mpirun -np 4 ./Fourier_Quad_Pipe /data/work/expo_gband.list
```

The Fortran executable accepts one positional exposure-list path; it has no
C++ phase CLI.

## Processing flow

The C++ driver can run five top-level phases:

```text
process_extcat -> process_init -> process_main -> process_rearr -> process_fd
```

The C++ and Fortran `process_main`/main programs implement nine prime-gated
numerical stages from preprocessing to catalog combination. See the dedicated
guides for input, configuration, and output contracts.

## Containers and HPC

- [`cpp_docker`](cpp_docker/) supplies the C++ toolchain image and a PMI2
  Slurm/Apptainer runner.
- [`f77_docker`](f77_docker/) supplies the pinned GNU 4.8.5/MPICH toolchain and
  generic plus pilogin Slurm templates.

Images contain toolchains and libraries only. Source, catalogs, observation
data, and outputs remain on bind-mounted host storage.

## Documentation

| Document | Purpose |
|---|---|
| [C++ guide](CPP_GUIDE.md) / [中文](CPP_GUIDE_CN.md) | Build, CLI, stages, inputs, and outputs |
| [C++ parameters](CPP_PIPELINE_PARAMETERS.md) | CLI/default mapping, compile-time controls, and 44-column default catalog |
| [F77 guide](F77_GUIDE.md) / [中文](F77_GUIDE_CN.md) | Fortran configuration, build, run, and outputs |
| [External catalogs](gen_src_cat/README.md) | Standalone repartitioner and TAP downloader |

## License

Repository-authored code is distributed under the [MIT License](LICENSE).
Container dependencies retain their upstream licenses; consult each container
directory's third-party notices.
