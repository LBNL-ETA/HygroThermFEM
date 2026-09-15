# Validation tests

Tests in this folder assert VALIDATED NUMBERS -- values whose authority lies
outside the engine: published analytic solutions (Carslaw & Jaeger), report
tables (UCRL-ID-106550), normative standards (EN 15026:2007 Annex A), the
HAMSTAD WP2 benchmark package (its analytic case, its participants' profiles
and its band of acceptance; data in hamstad/), or values independently
computed for the same case. Their purpose is regression tripwire: when
coupled or conduction physics drifts, these fail HERE, in the engine's own
suite.

Checkpoint tolerances are set at roughly twice the engine deviation measured
at first capture, except where a standard sets the band itself.

Membership criterion: a test belongs here if its EXPECT values come from a
publication, a standard, or an independent computation -- not if it only
checks internal consistency (conservation, symmetry, sanity bounds).

Two kinds of test live here. Most run a solve and compare the result. A few
instead check the MATERIAL a benchmark prescribes, against the values the
publication prints for it -- EN15026_AnnexAProperties is the example. Those
run in microseconds and fail with a specific message, where a transcription
error in a property function would otherwise surface as a vague disagreement
at the end of a long coupled solve, or not at all when the engine cannot
meet the publication's acceptance criterion for an unrelated reason.
