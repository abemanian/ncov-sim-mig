# ncov-sim-mig

Individual-based simulations (R package `nosoi`) for testing whether the normalized
relative risk (nRR) of identical pathogen sequences can detect changes in travel between
regions, and whether it stays robust when the force of infection changes.

## Model

Four demes in two regions: A–B (Region 1) and X–Y (Region 2), joined only by the A–X edge.

```
  B ---- A ---- X ---- Y
```

Each day a host moves with probability 0.10. The destination weights are:

|        | → A  | → B  | → X  | → Y  |
|--------|------|------|------|------|
| **A:** | —    | 0.70 | 0.30 | 0    |
| **B:** | 1.00 | —    | 0    | 0    |
| **X:** | 0.30 | 0    | —    | 0.70 |
| **Y:** | 0    | 0    | 1.00 | —    |

Movement restrictions scale per-edge movement rates. Suppressed travel is lost, not
redirected. The schedule has five phases:

| Phase | Start | All edges | A–X bridge (`heterogeneous`) |
|-------|-------|-----------|------------------------------|
| burn-in   | 0   | 1.00 | 1.00 |
| shutdown  | 100 | 0.10 | 0.10 |
| initial   | 150 | 0.33 | 0.10 |
| moderate  | 200 | 0.66 | 0.33 |
| full      | 250 | 1.00 | 0.50 |

Three movement scenarios (`--scenario`), each run with and without a transmission
intervention, give six arms:

| Scenario | Figure label | Movement |
|---|---|---|
| `none` | No Shutdown | baseline throughout |
| `homogeneous` | Uniform Recovery | schedule on every edge |
| `heterogeneous` | Staggered Recovery | same, but A–X recovers a phase late and only to 0.5 |

**Intervention** (`--vaccine`): from t = 175, transmission in Region 1 is multiplied by
0.8. Movement is unchanged.

### Parameters

All values are set as constants in `scripts/two_cluster/two_cluster_sim.R`. Time is in
days, and every host is updated once per day.

**Infection and transmission**

| Parameter | Value | Meaning |
|---|---|---|
| `INCUB_MEAN`, `INCUB_SD` | 3, 0.5 days | Incubation period, gamma with shape 36 and rate 12, drawn once per host; no transmission before it ends |
| `PROB_EXIT` | 1/7 per day | Daily probability that a host stops being infectious, applied from infection; mean time infected is 7 days |
| `CONTACT_MEAN`, `CONTACT_K` | 4, 0.2 | Contacts per day, negative binomial with mean 4 and dispersion 0.2 (strongly overdispersed) |
| `INF_ALPHA` | 20 | Shape α of the Beta distribution for per-contact transmission probability, redrawn each day |
| `INF_BETA_0` → `INF_BETA_STABLE` | 60 → 260 | Shape β, rising linearly over the seeding period and constant after it; mean transmission probability falls from 0.25 to 0.071 |
| `TIME_SEED` | 50 days | End of the seeding period |

**Intervention** (`--vaccine` arms only)

| Parameter | Value | Meaning |
|---|---|---|
| `TIME_VACC` | day 175 | Start of the intervention |
| `VE` | 0.2 | Relative reduction for transmission probability of hosts in Region 1 (A, B) after the intervention starts |

**Movement**

| Parameter | Value | Meaning |
|---|---|---|
| `PROB_MOVE` | 0.10 per day | Baseline probability that a host moves |
| `BASE_WEIGHTS` | table above | Destination probabilities given a move |
| `PHASE_STARTS` | 0, 100, 150, 200, 250 | Phase start days |
| `PHASE_ALL` | 1, 0.10, 0.33, 0.66, 1 | Multiplier on every edge's rate, per phase |
| `PHASE_BRIDGE_HETERO` | 1, 0.10, 0.10, 0.33, 0.50 | Multiplier on the A–X rate under `heterogeneous` |

**Run control**

| Parameter | Value | Meaning |
|---|---|---|
| `INIT_INDIVIDUALS`, `INIT_DEME` | 5, A | Number and deme of infected hosts at day 0 |
| `MAX_TIME` | 350 days | Run length |
| `MAX_OUTBREAK` | 10⁷ | Cap on total infections (never reached) |
| `MIN_INFECTED` | 200 | Minimum total infections for a run to be accepted |
| `MIN_INFECTED_BURN_IN` | 100 | Minimum infections before day 100 for a run to be accepted |
| `MAX_TRIES` | 50 | Attempts per seed before giving up |

Acceptance depends only on the burn-in, so all six arms accept the same attempt for a
given seed and stay paired.

`figs/pdf/movement_network.pdf` and `figs/pdf/movement_design.pdf` draw the network and
the schedules.

## Pipeline (per seed)

1. **Simulate** (`two_cluster_sim.R --skip_tree`): writes host metadata (`meta.tsv`).
2. **Subsample** 50,000 tips, split evenly across 10-day infection-time bins
   (`subsample_fasta.py --ids-out`).
3. **Build the tree** of the sampled tips from the transmission chain
   (`build_sampled_tree.R`). `validate_sampled_tree.R` checks it against nosoi's own tree.
4. **Evolve sequences** with phastSim: Wuhan-Hu-1 reference, rate 2 × 10⁻⁶ per site per
   day, seed 17.
5. **Find identical pairs** with `pairsnp -d 0`.
6. **Compute RR and nRR** (`two_cluster_RR.R`) per phase and in right-sided 40-day
   sliding windows ([t − 40, t], every 5 days):

   ```
   RR(X, Y)  = n_XY · N / (n_X · n_Y)
   nRR(X, Y) = RR(X, Y) / mean(RR(X, X), RR(Y, Y))
   ```

   Here `n_XY` is the number of identical pairs with one sequence in X and one in Y. nRR
   corrects for the regions' relative outbreak sizes, which shift raw RR even when
   movement is unchanged.
7. **Compute R(t) by region** (`two_cluster_reproduction_region.R`) from the full
   metadata.

## Getting Started

```bash
conda env create -f ncov-sim-env.yml
conda activate ncov-sim
Rscript scripts/common/install_packages.R
```

Run all six arms (about a day for 20 seeds per arm at 50,000 tips), then the combined
figures:

```bash
bash scripts/two_cluster/two_cluster_batch.sh all 20 50000 \
  > /dev/null 2>&1 & disown
```

Logs go to `results/two_cluster/logs/<date>/`, with overall progress in `driver.log`. To
run one arm, use `two_cluster_batch.sh [SCENARIO] [N_SEEDS] [SUBSAMPLE_N] [INTERVENTION]`.
For example, Staggered Recovery with the intervention:

```bash
bash scripts/two_cluster/two_cluster_batch.sh heterogeneous 20 50000 1
```

Run arms one at a time, because they share intermediate files. To redraw figures without
rerunning:

```bash
Rscript scripts/two_cluster/two_cluster_combined_plot.R   # figs/pdf/*_combined.pdf
Rscript scripts/two_cluster/two_cluster_network_plot.R    # figs/pdf/movement_*.pdf
```

The combined figures are laid out with movement scenarios as rows and intervention as
columns:

| File | Shows |
|---|---|
| `rt_region_combined` | R(t) by region |
| `rr_by_pair_combined` | sliding-window RR, one line per deme pair per seed |
| `nrr_by_pair_combined` | sliding-window nRR, one line per deme pair per seed |
| `recovery_by_phase_combined` | each phase's nRR relative to the [60, 100] pre-shutdown window |
| `recovery_by_phase_rr_combined` | the same for RR |

## Repository Layout

Run every script from the repository root.

```
scripts/
  two_cluster/   # sim, RR, R(t), batch, per-arm and combined plots, movement figures
  common/        # calc_rr.R, subsample_fasta.py, build/validate_sampled_tree.R,
                 # install_packages.R
  basic/         # single-population baseline
data/reference/  # wuhan.fasta
figs/pdf/        # combined and movement figures (tracked)
figs/jpg/        # per-arm diagnostic figures (gitignored)
results/         # per-arm RR / R(t) tables and logs (gitignored)
```
