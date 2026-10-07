# HygroThermFEM

Finite-element engine for two-dimensional heat and moisture transport. It solves steady-state
and transient problems on quadrilateral meshes, with thermal and moisture fields coupled through
temperature-dependent material properties, and is the solver behind THERM's hygrothermal runs.

The library ships **no mesher and no file format**. It takes a mesh, materials and boundary
conditions through its API and returns nodal results; loading THMZ archives and producing meshes
belong to the projects layered on top of it.

## Building

The project builds with CMake presets. `default-*` fetches every dependency from its declared
remote; `local-*` prefers a sibling working copy (`../<Name>`) of each dependency when one is
checked out, falling back to the remote per dependency, which is what you want while developing
across the dependency graph.

```
cmake --preset local-release
cmake --build --preset local-release
ctest --preset local-release
```

Swap `local-release` for `local-debug`, `default-release` or `default-debug` as needed. The build
directory is `build/<preset name>`.

`CMakeUserPresets.json` adds personal toolchain presets (gcc, clang, specific compiler versions,
Visual Studio) that inherit the same `local` behaviour.

### Dependencies

Fetched via `FetchContent`, pinned in `CMakeLists.txt`:

| Dependency | Pin | Purpose |
|---|---|---|
| Eigen | `5.0.1` | sparse matrix storage and linear solves |
| KeffCavity | `Version_1.1.5` | frame-cavity effective conductivity (also pins Windows-CalcEngine) |
| LBNLCPPCommon | `v0.18` | `ExpectedExt`/`OptionalExt` and the `lbnl::` algorithm helpers |
| GoogleTest | `v1.17.0` | tests only |
| FileParse | `Version_1.1.5` | tests only; loads the error-estimator regression fixtures |

### Build options

| Option | Default | Effect |
|---|---|---|
| `BUILD_HTFEM_TESTING` | `ON` standalone, forced `OFF` as a subproject | Builds `HygroThermFEM_tests` and `HygroThermFEM_sweep` |
| `HTFEM_ENABLE_LTO` | `ON` | Interprocedural optimization for Release; needed to inline `SquareMatrix` accessors into the assembly loop |
| `BUILD_FOR_CODE_COVERAGE` | `OFF` | Coverage instrumentation. GCC and Clang only; MSVC warns and builds uninstrumented |
| `BUILD_FOR_PROFILING` | `OFF` | Forces `-O0` for profile-friendly builds |

## Using it

Link the target; its public headers and Eigen's include path come with it.

```cmake
target_link_libraries(my_target PRIVATE HygroThermFEM)
```

```cpp
#include <HygroThermFEM/HygroThermFEM.hxx>
```

`MultiDomain` is the entry point: register materials, create elements over a node pool, attach
boundary conditions, then call `steadyState()` or `transient()`. Both return a `Solution` holding
nodal temperature, humidity, water/liquid/vapour/ice content and per-node fluxes.

Header layout:

- `include/HygroThermFEM/` -- curated facades: `Engine.hxx`, `Results.hxx`, `Config.hxx`,
  `Boundary.hxx`, `Materials.hxx`, `Geometry.hxx`, and `HygroThermFEM.hxx` pulling in all of them.
- `include/HygroThermFEM2D.hxx` -- back-compat umbrella exposing the full historical surface.
  Existing consumers use this; new code should prefer the facades.

## Testing

`HygroThermFEM_tests` is one GoogleTest binary, registered with CTest as a single entry. Run it
directly for filtering:

```
./build/local-release/Release/HygroThermFEM_tests --gtest_filter=SteadyState*
```

Tests live under `tst/units/`, grouped by `common/`, `steady_state/`, `transient/`, `enclosure/`,
`errorest/`, `therm_samples/` (regressions against THERM sample models) and `validation/`
(analytical and published-benchmark comparisons). `tst/sweep/` builds a separate
`HygroThermFEM_sweep` executable for solver parameter studies; it is not part of the suite.

## Layout

- `src/` -- implementation, grouped by component: `boundary/`, `common/`, `domain/`, `elements/`,
  `errorest/`, `materials/`, `math/`, `mesh/`, `simulation/`
- `include/` -- the public headers described above
- `tst/` -- tests and helpers
- `doc/` -- the numerical documentation

## Documentation

- [Discretization of governing equations](doc/Discretization%20Of%20Governing%20Equations.md)
- [Moisture governing equations](doc/Moisture%20Governing%20Equations.md)
- [Iteration methods](doc/IterationMethods.md)

## License

See [LICENSE](LICENSE) for details.
