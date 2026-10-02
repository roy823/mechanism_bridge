# FP-JCTC-1 (frozen 2026-10-02)

Frozen protocol for the JCTC campaigns: every SearchProtocol field except
random_seed, which is passed per run. FP-JCTC-1.json is byte-identical to the
last FP-JCTC-1-draft.json that passed the confirmation pilot. Each decision
follows a rule written before its test; numbers, job IDs and commits are in
docs/10-02 JCTC协议冻结前测试记录.md of the outer workspace.

| Field | Value | Evidence (UF HiPerGator, 2026-10-02) |
|---|---|---|
| ts_optimizer | dimer+sella | P1, 100 historical seeds: index-one TS 27% -> 44% at 0.96x median evaluations per TS |
| dimer_fmax, sella_steps | 1.0 eV/A, 400 | P1c, 120 seeds >= 25 atoms: accepted connections 1.7% -> 19.2%; 120 smaller seeds 47.5% -> 59.2% at 0.63x median evaluations |
| ts_steps | 160 | unchanged; a 400-step Dimer exhausts the attempt cap (P1c) |
| connection_protocol | irc | P1b: accepted 30% -> 37%; P2: endpoint graphs 43/43 identical; P2b: 7/7 DFT IRC endpoints at 1.59x endpoint evaluations (limit 2x) |
| fmax | 0.005 eV/A | P3: 0.003 versus 0.005 gives 50/50 identical TSs and endpoints |
| min_barrier_eV | 1e-3 | rigid-motion energy noise 1.35e-6 eV; rule max(1 meV, 5 sigma) |
| evaluations_per_attempt | 1900 | P90 of accepted attempts 1855 (P1b); kept under the P1c TS stage, where large-system accepted attempts reach a P90 of 1807 |
| encounter_policy | matched_controls | all four strategies get the same rigid encounter search |
| seed_features | arrow_features_v1 | pilot 3: lone-pair terms in 63% of arrows attempts, generalized links 100%, cyclic chains 34%, encounter alignment 38% |
| seed_fit_max_nfev | 3000 | the cap of 200 failed 20 of 260 arrows seed fits and no bond_edits fit; all 20 converge (nfev <= 998, median 3.9 s); converged fits unchanged |
| seed_feature_ablation | none | Fig. 5c arms change only this field |
| frontier_connectivity | species | pilot 3: every formaldehyde-dimer edge started from a re-oriented encounter complex of the root species, so physical connectivity never expanded the products |
| proposal_resonance_forms | 8 | product connectivity proposed for Coley 19/100 -> 86/100; T1x 13/204 unchanged |
| initial_steps | 1500 | Coley encounter roots stopped at 0.0099 eV/A after 250 BFGS steps |
| other fields | as in network_growth_v4 | 16000 evaluations, 24 attempts, batch-32 finite-difference Hessian |

Confirmation pilot v4 (UF job 44482971, f9548d6, 8 systems x 4 strategies x 2 seeds): 64/64 runs
completed with no seed-fit failure; symbolic runs on the Coley starts found 3-8 edges each
(pilot 3: almost none) and formaldehyde-dimer runs now have root-connected species pairs.

Run settings outside the protocol: AIMNet2-rxn on CPU, one thread and 2 GB per
run (P5), seeds 17 29 43 59 71 83 97 113 127 139 (P4 required more than five).

FP-JCTC-1-arrows-legacy.json differs only in seed_features='legacy'. It runs the
arrows strategy with the historical seed (0.7 source/sink coupling only) on the
Fig. 3/5a systems as the A_cur arm of the E-doc design; its seed RNG equals the
FP-JCTC-1 arrows arm, so the two are paired.

FP-JCTC-1-ablation-<arm>.json (Fig. 5c) each differ from FP-JCTC-1 only in
seed_feature_ablation; they run the arrows strategy on configs/campaigns/fig5c_ids.json.
