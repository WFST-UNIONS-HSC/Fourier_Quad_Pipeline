# cpp_Lite

Reduced C++17 Fourier_Quad pipeline for Gaia astrometry, no super-flat,
per-chip DQ masks, external sources, enabled deblending, local-polynomial PSF,
and no PCA. Alternate Standard branches are absent from this tree.

```bash
make -j4
./Fourier_Quad_Pipe --help
mpirun -np 4 ./Fourier_Quad_Pipe \
  --run-init false --run-main true --run-rearr false --run-fd false \
  --expo-list /data/work/expo_gband.list
```

Pass `STACK_PREFIX` and, when needed, `EIGEN_INCLUDE` for nonstandard library
locations. Current Make targets are `all`, `clean`, and
`test-general-infrastructure`.

Lite defaults to initialization and main enabled, with rearrangement and FD
disabled. Its removed branches cannot be enabled by adding constants.

Stage 7 writes 24 fields and Stage 9 appends exposure chi-square. See the
[C++ guide](../CPP_GUIDE.md) and
[parameter reference](../CPP_PIPELINE_PARAMETERS.md).
