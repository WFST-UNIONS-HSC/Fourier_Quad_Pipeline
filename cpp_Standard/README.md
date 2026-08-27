# cpp_Standard

Full C++17 Fourier_Quad pipeline with optional flat, mask, identity-astrometry,
external/hybrid PSF, and PCA/multi-scale branches.

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

Standard defaults to initialization, main, rearrangement, and FD enabled.
Review phase switches and paths before running. CLI overrides workflow paths;
scientific branches and thresholds in `config/LensingConfig.hpp` require
rebuilding.

Stage 7 writes 24 fields and Stage 9 appends exposure chi-square. See the
[C++ guide](../CPP_GUIDE.md) and
[parameter reference](../CPP_PIPELINE_PARAMETERS.md).
