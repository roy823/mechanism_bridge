# C2. Methods for the frozen protocol FP-JCTC-1 (JCTC manuscript draft)

- 版本：2026-10-02 起草，对应本地分支 `claude/jctc-code-fixes`。协议已于 2026-10-02 冻结：提交 9734566（本地标签 `protocol-FP-JCTC-1`），`configs/FP-JCTC-1.json` 的 SHA-256 为 56549ba7…a77fe。
- 与 `C_Methods初稿.md` 的关系：初稿基于 main@6715d82 的静态阅读，大部分"计划改动"现已实现，并经过冻结前测试（`docs/10-02 JCTC协议冻结前测试记录.md`）。本稿只写冻结协议的实际行为。初稿保留作历史记录。
- 标记：
  - [TBD-P1c]、[TBD-freeze]：已解决（P1c 采用 dimer_fmax 1.0、sella_steps 400；冻结信息见上）；
  - [VERIFY]：需要核对文献或库默认值。
- 正文用英文，便于直接进入稿件。

---

## 2. Methods

**Overview.** Each run starts from one user-supplied geometry, relaxed and curvature-checked on the machine-learned interatomic potential (MLIP) AIMNet2-rxn. At each selected minimum, a symbolic layer proposes two-electron actions (curly arrows). Each action is turned into a three-dimensional seed and an initial mode. A minimum-mode search hands over to a P-RFO refinement, the saddle point is tested with a finite-difference Hessian, and an IRC on the MLIP followed by a polish locates the two connected minima. Minima and transition states (TSs) are registered in a physical network, from which the next frontier node is chosen. Four seeding strategies share every step after the seed: `geometry`, `center_random`, `bond_edits` and `arrows`. They differ only in the information the seed may use. All settings form the frozen protocol FP-JCTC-1 (`configs/FP-JCTC-1.json`, SHA-256 56549ba7e86aaab94ec9d81aab1fd3a146655723f330b610e581e4a3f1fa77fe, code commit 9734566); every run manifest records this hash, the code commit and the runtime environment.

### 2.1 Symbolic action library

Published arrow-pushing records (SynEPD polar set) and a reactant-only analyst grammar of 14 named two-electron actions supply the proposals; both are described in ref. [prior work] and are unchanged here. Two additions affect FP-JCTC-1.

*Resonance-aware proposals.* Geometry perception assigns one Lewis structure, whereas templates and grammar rules are written for particular drawings. A nitrile ylide perceived as its allenyl–enolate form, for example, matches no dipole rule. Each node therefore also contributes proposals from up to eight resonance forms of its graph: RDKit `ResonanceMolSupplier` with unconstrained anions and cations and charge separation allowed, keeping only forms with the node's total charge. Forms keep the atom order and connectivity, so their arrows and bond edits apply to the same geometry. Proposals are deduplicated by (net bond edits, arrows) and record the form they came from. On the Coley [3+2] benchmark this raised the fraction of reference reactions whose product connectivity is proposed from 19% to 86% (Wilson 95% CI 78–92%). It made no difference on Transition1x, where 13 of 204 representable reactions are covered with or without resonance forms (Section 3.x).

*Given reaction centre (oracle) track.* For the engine-only comparison, `OracleReferenceLibrary` proposes exactly the reference reaction's net bond edits (Kekulé orders), and only at the root reactant. These runs know the answer and are reported separately.

### 2.2 Action scheduling

At a node, the next action minimizes, lexicographically: whether its predicted product graph is already known; the number of trials already spent on that product graph; the number of uses of the action; and its position in the diversified proposal list. Each action has three geometric variants *v* ∈ {0, 1, 2}. The variant sets the progress fraction (0.35/0.55/0.75), the asynchronous phase (0/−0.15/+0.15), the deformation norm (0.6/1.0/1.4 Å) and the angle priors together. `geometry` uses up to nine geometry seeds per node. The other strategies use symbolic actions only, one seed per (node, action, variant) visit.

### 2.3 Seed construction

All strategies share the following steps. The seed random-number generator is derived from the node graph, its rounded coordinates, the action's bond edits, the run seed and the variant, so matched strategies draw the same random numbers.

1. *Encounter orientation* (bimolecular nodes). A rigid search over 64 random fragment orientations places the forming pair 2.8 Å apart and penalizes clashes. Under FP-JCTC-1 (`encounter_policy='matched_controls'`), `geometry` and `center_random` receive the same search with a proposal-free contact pair, so all four strategies start from equally prepared encounter complexes.
2. *Targets.* For each edited pair, the target distance is the initial distance plus 1.1 Å for a broken bond, and *c*(*b*′)(*R*ᵢ + *R*ⱼ) otherwise, with *c* = 1.00, 0.88 and 0.80 for single, double and triple bonds. The normalized progress target is clip(*f*ᵥ + φᵥ sgn Δ*b*, 0.15, 0.90).
3. *Angle priors* derived from the edits alone: SN-type approach, H transfer and carbonyl addition.
4. *Least-squares deformation.* The residual collects edited distances, unchanged bonds, a Cartesian regularization, angle priors and clash penalties. It is minimized with SciPy `least_squares` (ftol = xtol = gtol = 10⁻⁵) for at most `seed_fit_max_nfev` = 3000 function evaluations. The historical cap of 200 failed 7.7% of arrow seeds in the pilot and no `bond_edits` seed. All of those fits converged within 1000 evaluations, and converged fits are unchanged by the larger cap.
5. *Normalization.* Rigid-body motion is projected out and the deformation is rescaled to the variant norm. The initial mode is a damped internal-coordinate tangent at the seed.

*Arrow-only features* (`seed_features='arrow_features_v1'`) are added for `arrows`; everything derivable from the net bond edits stays in the shared steps above:

- (a) **Lone-pair direction.** For an arrow that donates a lone pair of atom *d* into a new bond *d*–*a*, VSEPR lone-pair directions at *d* are computed from its neighbours. A hinge residual max(0, cos τᵥ − **l**·**ê**_{da})/0.3, with tolerance τᵥ = 10/20/30°, penalizes only misalignment of the closest lone pair beyond the tolerance. A donor with a single neighbour receives a cone-angle target instead (109.5/120/180° by steric number).
- (c) **Push–pull order.** Arrows are ranked by breadth-first order along the chain (an arrow precedes another when its head starts the next one). Each edited bond receives a progress shift λ(ρ − 0.5), where ρ is the normalized mean rank of the arrows touching it and λ ∈ {0, −0.3, +0.3} is paired with the variant as λ = Λ[(*v* + 1) mod 3]. Every pair of edited bonds touched by consecutive arrows is also coupled by the source–sink progress residual (weight 0.7).
- (f) **Encounter alignment.** The rigid encounter score gains 2 max(0, 0.8 − **l**·**ê**_{da})² per lone-pair-donating arrow.

*Matched control.* `bond_edits` receives the same progress-shift axis, with the bond order drawn from a random permutation on a separate random stream. Its proposals have the arrows removed before seeding, so it can never read them. In the pilot, the arrow features were active in the following fractions of `arrows` attempts: lone-pair terms 63%, generalized links 100%, cyclic chains 34%, encounter alignment 38%.

### 2.4 Saddle search and index check

The seed is refined by an ASE minimum-mode (Dimer) search started along the initial mode. Settings: separation 0.005 Å, maximum translation 0.1 Å, at most three rotations per step, rotational force thresholds 0.01/0.1 eV Å⁻¹. The search stops when the maximum per-atom force falls below `dimer_fmax` = 1.0 eV Å⁻¹ or after `ts_steps` = 160 steps. If the hand-off threshold is met, the point is refined with Sella P-RFO (order 1, Cartesian coordinates) to fmax = 0.005 eV Å⁻¹ within `sella_steps` = 400 steps.

- P1: on historical 9–15-atom seeds, the Dimer→Sella hand-off raised the index-one rate from 27% to 44% at equal cost.
- P1c: on 30–40-atom seeds, the 0.05 eV Å⁻¹ hand-off left most Dimer runs above the threshold, so Sella rarely ran. Handing over at 1.0 eV Å⁻¹ with up to 400 Sella steps raised the accepted-connection rate from 1.7% to 19.2% on 120 large seeds (index-one 1.7% to 40%) and from 47.5% to 59.2% on 120 smaller seeds, at 0.63× the median evaluations per attempt.

The Hessian is formed by central finite differences of forces (0.005 Å, batches of 32 geometries) and mass-weighted, and rigid modes are projected out. A point is an index-one saddle if its forces are converged and exactly one frequency lies below −30 cm⁻¹.

### 2.5 Endpoints

Both directions of the IRC are integrated on the MLIP with Sella's IRC (step 0.08 Å amu^½, at most 300 steps, fmax 0.05 eV Å⁻¹, inner fmax 0.02 eV Å⁻¹, at most 20 inner iterations). An inner-loop convergence failure stops that branch and is recorded as an unconverged IRC. Each branch end is polished with BFGS (maximum step 0.1 Å, fmax 0.005 eV Å⁻¹, at most 250 steps) and checked with the same Hessian test.

An attempt yields a validated connection if both endpoints are minima with converged forces and if E_TS − E_end ≥ 1 meV for both. The threshold is max(1 meV, 5σ), where σ = 1.35 µeV is the rigid-motion energy noise of AIMNet2-rxn measured on 60 structures × 20 random rigid motions.

Against seven DFT-verified events, IRC endpoints and the earlier mode-displacement descents both matched the DFT IRC endpoint graphs 7/7. The IRC costs a median 1.59× more endpoint evaluations, within the pre-registered 2× limit.

### 2.6 Event registration

Minima are perceived as Lewis graphs (RDKit `DetermineBonds`, system charge, charged fragments allowed); open-shell graphs are rejected. Two minima are the same node if their graphs are identical, their energies differ by ≤ 0.03 eV and their symmetry-aware RMSD is < 0.15 Å. A validated attempt with endpoint nodes B ≠ C creates the undirected edge B–C, whatever the source node and the predicted product were. The edge is `conformational` if B and C are resonance-equivalent and `chemical` otherwise. Its chemistry label (intermolecular heavy-atom bond, H transfer, intramolecular, …) is computed in the attempt's own atom frame.

### 2.7 Frontier selection

The root is relaxed with BFGS to 0.003 eV Å⁻¹ (at most 1500 steps) and curvature-repaired up to three times. The next node is drawn from the root-connected component. Unvisited species come first, then nodes with fewer visits.

Under FP-JCTC-1 (`frontier_connectivity='species'`), minima with the same graph count as connected. The reason is a bimolecular effect: every encounter is re-oriented, so the reactant side of a validated edge is usually a different orientation of the root complex than the root minimum itself. With physical connectivity, the pilot formaldehyde-dimer runs found the C–C coupling and H-transfer chemistry, but no edge touched the root minimum, and the products were never expanded.

### 2.8 Budget and failure accounting

A run stops at 16,000 MLIP evaluations, 24 attempts or an exhausted frontier, whichever comes first; one attempt may use at most 1900 evaluations. Under the adopted TS stage, accepted attempts on large systems reach a 90th percentile of 1807 evaluations. Evaluations are counted at the calculator, including the root relaxation, every Dimer, Sella, Hessian and IRC step, and failed attempts.

Failures are handled as follows:

- Recoverable failures end the attempt with a recorded status and the run continues. These are `ValueError`, `RuntimeError` and `ArithmeticError` raised by geometry perception, RDKit or the optimizers, including seed fits that reach the cap; CUDA and out-of-memory errors are excluded.
- Programming or infrastructure errors stop the run with `fatal_error` and are never silently skipped.

Runs use one CPU thread (P5: highest throughput per core).

### 2.9 Benchmarks

| Set | Source | Admitted | Representable on AIMNet2-rxn | Search starts | Role |
|---|---|---:|---:|---:|---|
| Transition1x test | official split, 287 reactions | 244 | 204 | 58 distinct reactants | main unimolecular set |
| Coley/Stuyver [3+2] | 1242 eligible CHNO, 17 strata, proportional sample | 100 | 73 | 100 | bimolecular set |
| RGD1 | 122,487 deduplicated pairs, 4 strata (b1f1, b2f2, b3f3+, other) | 100 | 70 | 70 | seen-domain positive control |

- *Admission.* Graphs of both endpoints are perceived from geometry; reactions whose connectivity does not change are rejected. For RGD1, the mapped CSV must also agree with the HDF5 reactant and product geometries.
- *Representability.* Each reference TS is re-optimized with Sella P-RFO on the MLIP (200 steps, screen-only budget of 4000 evaluations), checked for index one, and followed by the MLIP IRC. A reaction is representable if the IRC endpoints reproduce the reference connectivities under a reactant automorphism. Only representable reactions enter the main denominator. Coley author TS atoms are mapped onto our atom order through an element-labelled isomorphism of the TS connectivity: at a [3+2] TS every reactant bond survives. The mapping with the shortest forming bonds is chosen.
- *Shared reactants.* The 204 representable Transition1x reactions come from 58 reactant geometries. Blind searches therefore run once per reactant and are scored against every reaction of that reactant. Intervals use a cluster bootstrap over reactants.
- *Training overlap* (fixed before any scoring). Each reaction is labelled seen-reaction, seen-reactant, seen-formula or unseen against RGD1, which is part of AIMNet2-rxn's training data, using stereo-free InChIKeys. Transition1x: 19/120/65/0 of the 204; Coley: 0/1/2/97 of 100. Main conclusions use only the unseen and seen-formula layers. These labels are lower bounds, because the remaining AIMNet2-rxn training data are not public.

### 2.10 Recovery criteria

- *S1 (mapped).* An accepted MLIP edge whose two endpoint connectivities equal the reference reactant and product connectivities under one automorphism of the reactant graph.
- *S1 (unmapped).* The unordered pair of stereo-free endpoint SMILES equals the reference pair.
- *Evaluations to hit* counts initialization plus every attempt up to and including the first hit.

### 2.11 External baseline

Single-ended GSM (pyGSM) runs on the same MLIP through a budgeted ASE calculator; every evaluation, including aborted GSM runs, is billed against the same 16,000-evaluation budget. Driving coordinates come from one of three arms:

- the symbolic proposals' net bond edits (`pygsm_edits`);
- seeded random b2f2 rearrangements with only a connectivity-degree ceiling and a 4 Å formation cut-off (`pygsm_b2f2`);
- the reference edits (`pygsm_oracle`, which knows the answer).

GSM TS nodes go through the same Sella P-RFO, index check and IRC endpoint chain as our runs.

### 2.12 DFT verification

Sampled MLIP edges are verified at ωB97X/6-31G(d) (PySCF 2.10; grid level 3, SCF convergence 10⁻¹⁰ Eh, no density fitting). The workflow is: Sella P-RFO from the MLIP TS (fmax 0.01 eV Å⁻¹, at most 100 steps), an analytic Hessian, Sella IRC in both directions (at most 160 steps), BFGS polish, and endpoint Hessians and graphs. An edge is verified if the DFT TS has one imaginary frequency and both IRC endpoints are distinct minima with non-negative barriers.

The GPU route uses GPU4PySCF with the same functional, basis, grid and SCF settings. In Pdft-0, GPU and CPU energies agreed to 3 × 10⁻¹⁰ Eh (≤ 1.5 × 10⁻⁷ kcal mol⁻¹), and the GPU was 8–9× faster than 8 CPU threads for 10-atom TSs and 116× faster for a 40-atom encounter complex (SCF, gradient and analytic Hessian). A non-converged DFT calculation counts as "not certified", never as evidence against the reaction.

The Fig. 6 sample was fixed before any DFT calculation: 80 deduplicated edges, 16 per strategy stratum, random seed 20261002.

### 2.13 Statistics

Every strategy runs with 10 seeds (17, 29, 43, 59, 71, 83, 97, 113, 127, 139). The pre-registered power check (P4) on the pilot required more than five seeds for arrows versus bond_edits.

- Fig. 3/5a: the efficiency metric is the number of distinct root-connected chemical species pairs found within 16,000 MLIP evaluations; minima with the same graph form one species, so a reaction from a re-oriented encounter complex of the root counts. Per-system seed means are compared across systems with paired Wilcoxon signed-rank tests (arrows vs bond_edits; each strategy vs geometry). The anytime curves are read on a 500-evaluation grid, and their bands are 95% percentile bootstraps over systems. Runs with another protocol (the arrows-legacy arm, the ablations) form separate arms.
- Fig. 5b: for every group of proposals that share net bond edits but differ in arrows, the arrow-set × TS-cluster table of accepted attempts is tested by permutation (chi-square statistic, 10,000 draws). Each arrow set is compared with bond_edits by exact McNemar tests of acceptance and of intended success, paired by sample and seed. Holm corrections apply within each family across groups.
- Fig. 5c: each of six seed-feature ablations of the arrows strategy (10 systems × 10 seeds) is compared with the full arrows strategy by paired Wilcoxon tests over systems, with a Holm correction over the six ablations.
- Fig. 4: recovery rates carry 95% cluster-bootstrap intervals over starts (10,000 resamples, seed 20261002). A sensitivity analysis leaves out starts whose MLIP-relaxed root does not keep the reference reactant connectivity; for these starts the reference reaction can only be found after the network returns to the reference reactant.
- Fig. 6: confirmation rates carry Wilson 95% intervals; DFT non-convergence counts as not certified.
- Attempt-level diagnostics: every attempt's outcome (new or duplicate connection, same-basin return, unconverged TS force, failed calculation, budget exhaustion), its evaluations, and whether a connection joins the source graph to the proposal's predicted graph are reported per strategy.
- All decision rules, thresholds and sampling seeds were written before the corresponding runs and are recorded with the code.

### 2.14 Software

Python 3.12; NumPy 1.26.4; SciPy 1.15.3; ASE 3.29.0 [VERIFY]; RDKit 2024.03.2; NetworkX 3.2.1; PyTorch 2.8.0 (CPU); AIMNet 0.2.0 (AIMNet2-rxn weights, SHA-256 in each manifest); Sella 2.6.0; PySCF 2.10.0; GPU4PySCF 1.8.1 with CuPy 13.6.0 (CUDA 12.8, NVIDIA B200); pyGSM (commit recorded in each baseline manifest). Computations ran on UF HiPerGator (AMD EPYC CPU nodes; B200 GPU nodes for DFT).
