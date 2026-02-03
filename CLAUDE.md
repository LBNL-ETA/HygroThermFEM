# HygroThermFEM Code Review and Improvement Plan

## Overview

Analysis of the HygroThermFEM codebase - a C++20 Finite Element Method engine for hygrothermal (coupled heat and moisture) analysis.

**Codebase Statistics:**
- 74 source files (~11,678 lines of C++ code)
- Dependencies: Eigen 5.0.1, KeffCavity 1.0.31, Google Test
- Build system: CMake with FetchContent

---

## Code Quality Issues

### 1. Variable Naming Conventions

**Problem:** Some single-letter variable names reduce readability.

**Note:** Standard FEM/linear algebra conventions should be preserved:
- `M` - Mass matrix (standard)
- `K` - Stiffness/conductance matrix (standard)
- `C` - Capacitance matrix (standard)
- `A` - Left-hand side matrix in Ax=B (standard)
- `B` - Right-hand side vector in Ax=B (standard)

**Variables to rename (non-standard single letters):**
```cpp
// Examples that should be renamed
double h;           // -> boundaryMatrix or heatTransferCoeff
double k;           // -> conductivity (when not stiffness matrix)
double a, b, c;     // -> extrapolationCoeffA, extrapolationCoeffB, etc.
auto cf;            // -> coefficient
auto bc;            // -> boundaryCondition
```

**Affected files:**
- `src/Domain.cxx` - `h` variable for boundary matrix
- `src/Element2D.cxx` - extrapolation coefficients `a`, `b`, `c`
- `src/BoundaryCondition2D.cxx` - temperature variable `T`

**Recommendation:** Rename non-standard variables to be at least 3 characters with descriptive names.

---

### 2. Code Duplication in Boundary Condition Classes

**Problem:** Significant code repetition across BC implementations.

**Duplicated patterns:**
- `ThermalBCConvection` / `MoistureBCConvection`
- `ThermalBCRadiation` / similar moisture classes
- `ThermalBCFlux` / `MoistureBCFlux`

**Recommendation:** Extract common BC logic into a template base class:
```cpp
template<typename DomainType, typename CoeffType>
class ConvectionBCBase {
protected:
    void applyConvection(DomainType& domain, CoeffType coefficient);
};
```

---

### 3. Const-Correctness Issues

**Problem:** Many getter methods and operators are not marked `const`.

**Examples found:**
```cpp
// Current
double getValue() { return value; }  // Should be const
Material& operator[](size_t idx);    // Should have const overload

// Recommended
double getValue() const { return value; }
const Material& operator[](size_t idx) const;
Material& operator[](size_t idx);
```

**Affected areas:**
- MaterialPool accessors
- Element property getters
- Domain query methods

---

### 4. MaterialPool Singleton Coupling

**Problem:** Global singleton `MaterialPool::Instance()` creates tight coupling and makes testing difficult.

**Current usage:**
```cpp
// Scattered throughout codebase
auto& material = MaterialPool::Instance().getMaterial(id);
```

**Recommendation:** Inject MaterialPool as dependency:
```cpp
class Domain {
public:
    explicit Domain(MaterialPool& materials);
private:
    MaterialPool& materialPool;
};
```

---

## Architectural Improvements

### 1. Boundary Condition Factory Pattern

**Problem:** BC creation is scattered with switch statements or if-else chains.

**Recommendation:** Implement factory pattern:
```cpp
class BCFactory {
public:
    static std::unique_ptr<IBoundaryCondition> create(
        BCType type,
        const BCParameters& params
    );

    // Registration for extensibility
    static void registerBC(BCType type, CreatorFunc creator);
};
```

---

### 2. Builder Pattern for Complex Objects

**Problem:** MultiDomain and GlazingSystem construction involves many parameters.

**Recommendation:** Implement builder pattern:
```cpp
auto system = MultiDomainBuilder()
    .addThermalDomain(thermalConfig)
    .addMoistureDomain(moistureConfig)
    .withCoupling(CouplingType::Full)
    .withMesh(mesh)
    .build();
```

---

### 3. Strategy Pattern for Solvers

**Problem:** Solver selection logic embedded in domain classes.

**Recommendation:** Extract solver strategies:
```cpp
class ISolverStrategy {
public:
    virtual void solve(SystemMatrix& matrix, Vector& rhs, Vector& solution) = 0;
};

class DirectSolver : public ISolverStrategy { /* LU decomposition */ };
class IterativeSolver : public ISolverStrategy { /* CG, GMRES */ };
```

---

### 4. Error Handling Improvements

**Problem:** Mix of exceptions and error codes, some errors silently ignored.

**Recommendation:**
- Use `std::expected<T, Error>` (C++23) or similar pattern
- Create domain-specific exception hierarchy
- Add error context for debugging

---

## Specific Refactoring Tasks

### High Priority ✓ COMPLETED

1. ✓ **Rename non-standard single-letter variables** - Completed
   - Renamed `h` → `boundaryHMatrix`, `k` → `conductivityValues`
   - Renamed `a`, `b`, `c` → `extrapolationCoeffA/B/C`
   - Renamed `T` → `nodeTemperature`, `surfaceTemp`, `radiationTemp`

2. ✓ **Add const to getters** - Completed
   - Added const to 8 methods: getLumpedMass, getMassMatrix, getSide1, getSide2,
     Psi, PsiDKsi, PsiDEta, DpDuMatrices

3. ✓ **Extract BC base classes** - Completed (3 phases)
   - Phase 1: Created `IRadiationBC` base for BlackBodyRadiationBC and LinearizedRadiationBC
   - Phase 2: Created `IEnvironmentCoefficients` base struct for 6 coefficient structs
   - Phase 3: Created `IConvectiveBCBase` for IConvectionBC and IMoistureBC

**Additional completed work:**
- Removed dead code: checkSingularity(), degrees(), State operators, TimestepData setters
- Fixed bug: ThermalConductivityDry missing from isMissingAnyProperty()
- Added unit tests: Domain_ClearModelAndGravity (5 tests), MaterialMissingPropertiesMessage (13 tests)

### Medium Priority (Next)

4. **Inject MaterialPool** - Remove singleton dependency
5. **Implement BC Factory** - Centralize BC creation
6. **Add missing unit tests** - Coverage for edge cases

### Lower Priority

7. **Builder pattern for domains** - Improve construction API
8. **Solver strategy extraction** - More flexible solver selection
9. **Documentation improvements** - Add Doxygen comments

---

## Testing Improvements

**Current state:** Google Test framework in place with basic coverage.

**Recommendations:**
1. Add parameterized tests for element types
2. Add integration tests for coupled analysis
3. Add benchmark tests for performance regression
4. Mock MaterialPool for isolated unit tests

---

## Performance Considerations

1. **Matrix assembly** - Consider parallel assembly for large meshes
2. **Memory allocation** - Pre-allocate vectors where sizes are known
3. **Eigen expressions** - Ensure lazy evaluation is preserved (avoid `.eval()` where not needed)

---

## Summary

The HygroThermFEM codebase has a solid foundation with proper C++20 usage and clean separation of thermal/moisture domains. The main areas for improvement are:

1. **Readability**: Variable naming (non-standard abbreviations) and const-correctness
2. **Maintainability**: Reduce BC code duplication
3. **Testability**: Decouple MaterialPool singleton
4. **Extensibility**: Factory and strategy patterns

These changes would improve code quality without requiring architectural rewrites.

---

## Related Libraries

### LBNLCPPCommon
**Location:** `D:\Programming\GitHub\LBNLCPPCommon`

Header-only C++20 utility library from LBNL. Provides generic algorithms and functional programming utilities.

**Components:**
| Header | Purpose |
|--------|---------|
| `algorithm.hxx` | `find_element()`, `filter()`, `transform_to_vector()`, `zip()`, `split()`, `flatten()` |
| `optional.hxx` | Extended optional with monadic ops: `and_then()`, `or_else()`, pipe operators (`\|`, `\|\|`) |
| `expected.hxx` | Result type for error handling without exceptions |
| `map_utils.hxx` | `map_lookup_by_key()`, `map_lookup_by_value()`, `map_keys()`, `map_values()` |
| `enum_index_mapper.hxx` | Bidirectional enum-to-index mapping |
| `memoize.hxx` | Thread-safe memoization (LazyEvaluator) |

**Namespace:** `lbnl::`

**Usage:** Currently used by LibraryFEMTherm for database lookups and optional chaining. HygroThermFEM could adopt for similar patterns.

**Integration via CMake:**
```cmake
FetchContent_Declare(
    LBNLCPPCommon
    GIT_REPOSITORY https://github.com/LBNL-ETA/LBNLCPPCommon.git
    GIT_TAG v0.14
)
```

---

### LibraryFEMTherm
**Location:** `D:\Programming\GitHub\LibraryFEMTherm`

LBNL C++20 library for thermal simulation data management. Depends on LBNLCPPCommon.

**Provides:**
- Materials database (solid, cavity, radiation enclosure)
- Gases database (pure gases, mixtures)
- Boundary conditions (steady-state and transient)
- THMZ file format (compressed XML)
- Mesh data structures

**Key namespaces:** `MaterialsLibrary`, `GasesLibrary`, `BCSteadyStateLibrary`, `BCInputFileLibrary`, `ThermFile`

**Dependencies:** LBNLCPPCommon, FileParse, miniz

**Potential integration:** Data layer for material/BC management, THMZ file I/O.

---

## Consumer: Therm Application

**Location:** `D:\Programming\GitHub\Windows-Tools\Products\Therm\source`

Windows application using HygroThermFEM as solver engine. Interface changes require coordinated updates.

### Files with HygroThermFEM Includes

| File | Purpose |
|------|---------|
| `BCInputFileRecords.h` | BC coefficient storage |
| `HygroThermFEM/HygroThermFEMGeometryBuilder.cpp/h` | CON file to MultiDomain adapter |
| `HygroThermFEM/HygroThermFEMBCBuilder.cpp/h` | BC creation adapter |
| `HygroThermFEM/HygroThermFEMResults.cpp/h` | Solution extraction |
| `Simulate.cpp/h` | Transient loop driver |
| `ThermDoc.cpp/h` | Main document class |
| `Dialogs/ProgressDialogHygroThermFEM.h` | Progress observer |
| `Dialogs/MissingMaterialPropertiesDlg.cpp` | Validation UI |

### Public API Used by Therm

#### MultiDomain Class Methods
```cpp
void clearModel()
MaterialsErrorCheckVector checkForMaterialsValidity(SimulationType type)
void performMoistureSimulation(bool enable)
void performThermalSimulation(bool enable)
Solution transient(vector<double>& temperatures, vector<double>& humidities, double dTime, size_t timestepIndex)
void subscribeMoisture(TimestepObserver* observer)
void unsubscribeMoisture(TimestepObserver* observer)
void subscribeThermal(TimestepObserver* observer)
void unsubscribeThermal(TimestepObserver* observer)
```

#### Boundary Condition Creation Methods
```cpp
void createBC_TARPHc(INode2D& node1, INode2D& node2, vector<TARPCoefficients>& coeffs, double surfaceTilt)
void createBC_ASHRAEInsideHc(INode2D& node1, INode2D& node2, vector<ASHRAEInsideCoefficients>& coeffs, double length, double tilt)
void createBC_ASHRAEOutsideHc(INode2D& node1, INode2D& node2, vector<ASHRAEOutsideCoefficients>& coeffs)
void createBC_FixedHc(INode2D& node1, INode2D& node2, vector<FixedBCHCCoefficients>& coeffs)
void createBC_YazdanianKlemsHc(INode2D& node1, INode2D& node2, vector<YazdanianKlemsCoefficients>& coeffs)
void createBC_KimuraHc(INode2D& node1, INode2D& node2, vector<KimuraCoefficients>& coeffs)
void createBC_LinearizedRadiation(INode2D& node1, INode2D& node2, vector<LinearizedRadiationBCCoefficients>& coeffs)
void createBC_BlackBodyRadiation(INode2D& node1, INode2D& node2, vector<BlackBodyRadiationBCCoefficients>& coeffs)
void createBC_FixedHeatFlux(INode2D& node1, INode2D& node2, vector<double>& flux)
void createBC_FixedTemperature(INode2D& node1, INode2D& node2, vector<double>& temperature)
void createBC_FixedTemperatureAndHumidity(INode2D& node1, INode2D& node2, vector<TemperatureAndHumidity>& values)
```

#### NodePool Singleton
```cpp
static NodePool& Instance()
INode2D& createNode(size_t nodeID, double x, double y, State initialState)
vector<double> properties(Variable var)  // Variable::temperature or Variable::humidity
INode2D& getNode(size_t nodeID)
```

#### MaterialPool Singleton
```cpp
static MaterialPool& Instance()
IMaterial& createSolidMaterial(string_view name)
IMaterial& createGas(string_view name, CavityStandard standard, GasProperties props)
```

#### SimulationProperties Singleton
```cpp
static SimulationProperties& Instance()
void setCalculationParameters(bool excludeLiquid, bool excludeEvapHeat, bool excludeCapillary, bool excludeVaporDiff, bool excludeTempMoistureDep)
void setIterationParameters(double relaxation, double tolerance, size_t maxIterations)
```

#### Data Structures
```cpp
struct Solution {
    vector<double> temperature;
    vector<double> humidity;
    vector<double> waterContent;
    vector<NodeFlux> heatFlux;
    vector<NodeFlux> waterFlux;
    double dTime;
};

struct State {
    double temperature;
    double humidity;
};

struct NodeFlux {
    double x;
    double y;
};

struct TemperatureAndHumidity {
    double temperature;
    double humidity;
};
```

#### Enumerations
```cpp
enum class SimulationType { Thermal, Moisture, HygroThermal };
enum class Variable { temperature, humidity };
enum class CavityStandard { ISO15099 };
enum class WindDirection { Windward, Leeward };
```

#### Coefficient Structures (all used in vectors for timestep data)
- `TARPCoefficients`
- `ASHRAEInsideCoefficients`
- `ASHRAEOutsideCoefficients`
- `FixedBCHCCoefficients`
- `YazdanianKlemsCoefficients`
- `KimuraCoefficients`
- `LinearizedRadiationBCCoefficients`
- `BlackBodyRadiationBCCoefficients`

### Interface Stability Rules

**BREAKING (requires Therm updates):**
- Any signature change to methods listed above
- Removing or renaming public methods
- Changing data structure field names/types
- Changing enum values
- Removing singleton Instance() patterns

**SAFE (no Therm impact):**
- Internal/private implementation changes
- Adding new public methods (additive)
- Adding const qualifiers to existing methods
- Performance optimizations
- New optional parameters with defaults
