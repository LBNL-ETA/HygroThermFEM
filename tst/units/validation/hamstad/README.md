# HAMSTAD WP2 benchmark data

Input and reference data of the HAMSTAD WP2 benchmark package, read at runtime by
`tst/units/validation/HAMSTAD_*.unit.cxx` through the `HAMSTAD_DATA_DIR` compile
definition. The benchmark descriptions themselves (geometry, material functions,
boundary and initial conditions) are transcribed into `tst/helper/HAMSTADBenchmarks.hxx`
with page references; only what the package distributes as data lives here.

Source: the HAMSTAD project archive as redistributed by SimQuality
(https://simquality.de/wp-content/uploads/2021/05/HAMSTAD-Project.zip), containing
Hagentoft, C-E. 2002, *HAMSTAD - Final report: Methodology of HAM-modeling*, Report
R-02:8, Chalmers University of Technology, together with the benchmark descriptions,
climate files and the participants' result spreadsheets. The benchmarks were published
for exactly this use: the assessment of heat, air and moisture simulation codes.

| File | Content | Origin in the archive |
| --- | --- | --- |
| `ClimateBench1.txt` | Benchmark 1 climate, hourly for one year: time (s), T_eq,e (C), T_eq,i (C), p_a,e (Pa), p_a,i (Pa) | `HB1/ClimateBench1.txt`, verbatim |
| `Bench1Year1_MA.csv`, `Bench1Year1_MB.csv` | Benchmark 1, first year: integrated moisture (kg/m2) in the load bearing layer (MA) and the insulation (MB), hourly, per participating code, with the average, standard deviation and the 99.9 % band of acceptance | `HB1/Bench1Year1.ods`, sheets `ResultsMA` / `ResultsMB`, verbatim |
| `Bench1Year5_MA.csv`, `Bench1Year5_MB.csv` | The same for the fifth year | `HB1/Bench1Year5.ods`, sheets `Results5MA` / `Results5MB`, verbatim |
| `Bench5_<code>.csv` | Benchmark 5 at 60 days: position (mm), water content (kg/m3) and relative humidity (-) as submitted by CTH, KUL, NRC, Technion, TUD and TUE | `HB5/Bench5.xls`, sheet `Data`, one participant's final columns each (TUD's base run, Technion's final submission, KUL's `w (KUL)` / `fi (KUL)`, TUE's fractional humidity column, CTH's `w(CTH)` / `fi(CTH)`, NRC) |

Benchmark 2's analytical solution is a 21-by-3 table and is transcribed directly into
its test. Benchmarks 3 (air transfer through the layer) and 4 (driving rain with a
saturation cap) need physics the engine does not carry and are not included.
