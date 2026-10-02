# FP-JCTC-1 (draft, 2026-10-02)

Candidate frozen protocol for the JCTC campaign; every SearchProtocol field
except random_seed (passed per run). Not yet frozen: the pilot campaign must
confirm the arrow features on real systems, and P2b (MLIP IRC versus DFT IRC
endpoints) decides connection_protocol.

| Field | Value | Evidence (UF runs, 2026-10-02) |
|---|---|---|
| min_barrier_eV | 1e-3 | rigid-motion noise sigma 1.4e-6 eV (aimnet2-rxn), 1.5e-6 eV (aimnet2-2025); rule max(1 meV, 5 sigma) |
| ts_optimizer | dimer+sella | P1, 100 replayed seeds: index-one TS 27% -> 44%, median evaluations per TS 0.96x; rule (2) met |
| connection_protocol | irc (pending P2b) | P1b: accepted connections 30% -> 37% at equal median cost; P2: endpoint graphs 43/43 identical where both accepted |
| evaluations_per_attempt | 1900 | P90 of accepted attempts under dimer+sella+irc = 1855, rounded up to 100, floor 1400 |
| encounter_policy | matched_controls | fair rigid-encounter search for the geometry/center_random controls |
| seed_features | arrow_features_v1 | arrow-only lone-pair, push-pull order and coupling terms; matched random-order control for bond_edits |
| other fields | as in v4 (network_growth_v4) | 16000 evaluations, 24 attempts, fmax 0.005 eV/A, batch-32 FD Hessian |
