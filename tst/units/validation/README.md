# Validation tests

Tests in this folder assert VALIDATED NUMBERS -- values whose authority lies
outside the engine: published analytic solutions (Carslaw & Jaeger), report
tables (UCRL-ID-106550), or the externally validated 1D reference solver
(hygrotherm1d, the hygrothermfem_python validation book, which itself passes
EN 15026:2007 Annex A). Their purpose is regression tripwire: when coupled or
conduction physics drifts, these fail HERE, in the engine's own suite,
without needing the Python side.

The full-profile comparisons, provenance and tolerance derivations live in
the validation book (hygrothermfem_python); each dump-producing test names
its book dataset. Checkpoint tolerances are set at roughly twice the engine
deviation measured at first capture.

Membership criterion: a test belongs here if its EXPECT values come from a
publication, a standard, or the reference solver -- not if it only checks
internal consistency (conservation, symmetry, sanity bounds).
