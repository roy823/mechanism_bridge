# Architecture

## Scientific flow

```text
Observed reactant + public symbolic proposals
 -> geometric seed and internal direction
 -> unbiased MLIP saddle search
 -> Hessian and two-sided descent
 -> actual endpoint registry and parallel TS edges
 -> reachable-frontier exploration
 -> selected DFT/IRC verification
 -> electronic interpretation and reviewed paired data
```

| Concern | Modules | Contract |
|---|---|---|
| Source data and graphs | adapters, chemistry, event_graph, io | Atom identity, explicit system and provenance |
| Symbols | symbolic_library, local_patterns | Occupied electron-source checks; charge/valence local SMARTS; preserve distinct source/sink alternatives |
| Seeds | search_seeds | Equal displacement norms; asynchronous distance/angle priors and seed-local internal-coordinate tangent; no reference TS or product geometry |
| Energy/forces | aimnet_backend, backends | Official model or explicit PySCF; separate energy scales |
| Local physics | physics | Projected Hessian and explicitly named mode-displacement descent |
| Exploration | reaction_network | Budgets, Dimer, observed minima, parallel edges and reachable frontier |
| DFT verification | verification | TS refinement, strict Sella IRC, endpoints; optional separate orbital annotation |
| Electronic interpretation | electronic, event_graph | IAO/IBO and candidate arrows; no independent-truth claim |
| Visualization | molecular_visuals | Saved coordinates/energies, 2D depictions, offline WebGL and PNG/PDF |
| Historical support | pairing, training_data, kinetics | Source matching, reactive training samples and explicit free-energy kinetics |

## Entrypoints

- `scripts/data`: downloading, importing and preparing initial states.
- `scripts/exploration`: individual runs, configured campaigns, audits, summaries and visualization.
- `scripts/qc`: DFT/IRC and orbital analysis.
- `scripts/baselines`: earlier RF and residual-update probes.
- `scripts/diagnostics`: data and execution diagnostics.
- `mechbridge` CLI: ingestion, pairing/splitting, XYZ export, IBO and local downhill checks.

Scripts order experiments; reusable scientific functions live in `src/mechbridge`. Former script paths have no compatibility aliases.

## Evidence and reproducibility

New runs store protocol, starting records, source ZIP and hashes, model and symbol-library hashes, package versions and every attempt's coordinates/trajectories. Budgets include failures and finite-difference checks.

Edges connect observed minima, independently of the node that proposed a seed. Root-connected edges, disconnected discoveries, chemical graph pairs and conformers have separate counts.

v1 is historical evidence. Its source is archived in `reports/repository_snapshots/pre_validation_v2.zip`. v2 changes seed normalization and adds a center-only control; current scripts do not reproduce v1's exact numbers. Raw data, weights, caches and newly generated reports are ignored by Git; already tracked historical artifacts remain preserved.

v3 replaces whole-reactant matching with local patterns and changes seed geometry/directions. Historical source ZIPs remain the evidence for earlier runs. Root initialization uses tighter force convergence and bounded negative-mode displacement to escape stationary torsional maxima; its evaluations are charged to the same run budget. Candidate paths still use the documented MLIP mode-displacement descent, not DFT IRC.
