# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repository contains epidemiological simulations using the `nosoi` package in R to model pathogen transmission across structured populations. The simulations generate transmission trees and phylogenetic data, which are then analyzed using sequence evolution tools (phastSim) and pairwise SNP analysis (pairsnp).

## Environment Setup

**Conda environment**: Use `ncov-sim` environment (Python 3.12 required)

```bash
conda env create -f ncov-sim-env.yml
conda activate ncov-sim
```

**R package installation**:
```bash
Rscript scripts/common/install_packages.R
```

Required packages include: nosoi, tidyverse, viridis, igraph, ggraph, ape, tidytree, argparse, treeio (Bioconductor), ggtree (Bioconductor)

## Repository Layout

Scripts live under `scripts/`, one subfolder per scenario plus a shared `common/`.
Simulation *outputs* live under `results/`, one subfolder per scenario. `data/` holds
only static inputs -- currently `data/reference/` with the Wuhan-Hu-1 reference genome.

**All scripts assume the repository root as the working directory.** Paths inside the
scripts (`results/<scenario>/...`, `source("scripts/common/calc_rr.R")`) are relative to the
root, so invoke them as `Rscript scripts/<scenario>/<script>.R`, never by `cd`-ing into
the scenario folder.

## Simulation Architecture

### Core Simulation Types

All simulations use the `nosoi` package to model individual-based transmission
dynamics. Apart from `basic_sim`, each places hosts in spatial demes connected by a
movement matrix -- three demes (A, O, B) in the threebody models, six in `two_cluster`.

1. **scripts/basic/basic_sim.R** - Single-population baseline simulation with no spatial structure
2. **scripts/shutdown_threebody/shutdown_threebody_sim.R** - Models intervention scenarios with time-varying contact rates and movement restrictions
3. **scripts/carrying_capacity/carrying_capacity_sim.R** - Includes density-dependent transmission via logistic carrying capacity, plus optional seasonal forcing, isolation, and immunization scenarios
4. **scripts/two_cluster/two_cluster_sim.R** - Six demes as two fully-connected
   triangles (A-B-C and X-Y-Z) joined by the single A-X edge, which is the only
   inter-cluster route. Asks whether an edge-specific travel policy is detectable in
   genetic data: a five-phase movement schedule (burn-in / shutdown / initial /
   moderate / full) is applied to every edge under `--scenario homogeneous`, while
   `heterogeneous` runs the same schedule but holds the A-X bridge a phase behind,
   never restoring it past half of baseline. `none` is the null model. Suppressed
   travel is lost rather than redirected, so within-cluster rates are unaffected by
   bridge policy.

   Two details worth knowing before editing this script: nosoi fixes
   `structure.matrix` for a whole run, so time-varying movement is achieved by
   patching `nosoi:::moveFunction.discrete` at load time (`patch_move_function()`,
   which fails loudly if the upstream signature changes). And acceptance in
   `simulate_until_size()` is judged on the burn-in period *only*, so all three
   scenarios accept the same attempts and stay paired by seed.

### Common Simulation Parameters

All simulations define:
- `p_Exit_fct`: Probability of host exit (death/recovery)
- `n_Contact_fct`: Number of contacts per time step (Poisson, or negative binomial where overdispersion matters, as in `two_cluster`)
- `p_Trans_fct`: Transmission probability (often beta-distributed, may depend on incubation period and time)
- `d_Incub`: Incubation period distribution (gamma)
- `p_Move_fct`: Movement probability between demes (structured models only)

Key simulation controls:
- `MAX_OUTBREAK`: Maximum number of infections
- `MAX_TIME`: Maximum simulation time steps
- `MIN_INFECTED`: Minimum outbreak size for acceptance
- `simulate_until_size()`: Wrapper that retries simulations until minimum size reached

### Simulation Outputs

Each simulation produces:
- `tree.nwk` or `*.nwk`: Newick format transmission tree
- `meta.tsv` or `*_meta.tsv`: Metadata including infection times, locations, and host IDs
- Optional plots (when `--draw` flag used)

## Analysis Pipeline

### Standard Workflow

1. **Run simulation** (e.g., `Rscript scripts/shutdown_threebody/shutdown_threebody_sim.R --seed 42`)
2. **Sequence evolution** via phastSim (molecular clock simulation on transmission tree)
3. **Subsampling** (optional) via `scripts/common/subsample_fasta.py`, which thins the
   phastSim FASTA to a target tip count before the O(n^2) pairsnp step. With `--meta`
   it splits the budget *evenly across infection-time bins* rather than uniformly over
   the file; uniform sampling tracks incidence and leaves the early epidemic with too
   few tips for a stable RR estimate. Sampling stays uniform within each bin, so deme
   composition inside a time window is unbiased. Used by `two_cluster`.
4. **Pairwise SNP analysis** via pairsnp (identify identical/similar sequences)
5. **Relative risk calculation** using `scripts/common/calc_rr.R` (geographic association analysis)

### Batch Execution

Shell scripts run parameter sweeps across multiple seeds:
- `scripts/shutdown_threebody/shutdown_threebody_batch.sh` - Shutdown scenario iterations (50 seeds)
- `scripts/carrying_capacity/carrying_capacity_batch.sh` - Carrying capacity scenario iterations
- `scripts/two_cluster/two_cluster_batch.sh [SCENARIO] [N_SEEDS] [SUBSAMPLE_N]` -
  Two-cluster sweep, defaults `homogeneous 25 8000`. Run it once per scenario; it
  clears only the accumulators for the scenario it was given, and every output file is
  prefixed with the scenario name so sequential sweeps do not overwrite each other.
  Ends by invoking `two_cluster_plot.R` for that scenario.

Pattern: Each script loops over seeds, runs simulation → phastSim → (subsample) → pairsnp → RR calculation → cleanup

### Command Line Arguments

Most R scripts accept:
- `--seed INTEGER`: Random seed (default: 66)
- `--draw`: Generate diagnostic plots (default: FALSE)
- `--scenario STRING`: Scenario variant, where the model defines them. The
  `two_cluster` scripts take `none`, `homogeneous` (default), or `heterogeneous`, and
  the value must be passed consistently to `_sim.R`, `_RR.R`, and `_plot.R` -- it
  selects the intervention *and* namespaces the output files.

Example:
```bash
Rscript scripts/shutdown_threebody/shutdown_threebody_sim.R --seed 123 --draw
Rscript scripts/two_cluster/two_cluster_sim.R --seed 123 --scenario heterogeneous --draw
```

## Key Functions and Utilities

### scripts/common/calc_rr.R

Calculates geographic relative risk (RR) from pairwise identical sequences:
- `calc_rr()`: Computes RR matrix for sequence pairs between demes
- `normalized_rr()`: Normalizes RR by diagonal entries (nRR metric)
- Supports time window filtering for temporal analysis

RR interpretation: RR(X,Y) > 1 indicates more genetic similarity between demes X and Y than expected by chance

`two_cluster_RR.R` runs this at two resolutions -- the 6x6 deme matrix and a collapsed
C1/C2 cluster matrix -- and on two time decompositions: one estimate per intervention
phase (pooling all pairs in the phase, which is far more stable than taking a median
over narrow windows) and a 20-step sliding window every 5 steps.

The `recovery_by_phase` plot does *not* use the burn-in phase as its baseline. It divides
each pair's phase nRR by that pair's sliding-window nRR at `tpoint = 90` (the window
[70, 110]), which is wide enough to stay stable for the sparse inter-cluster pairs. The
narrower alternatives were tested and rejected: a strictly pre-intervention window loses
too many inter-cluster pairs to emptiness, and a later-centred one is deflated by
shutdown-era transmission, which pushes every recovery ratio above 1. Burn-in itself is
dropped from the figure -- the seeding period is intrinsically dense in identical
sequences and not comparable to the later phases -- so `*_recovery_by_phase.tsv` now holds
four phases, not five.

That figure is drawn on a pseudo-log axis (`PSEUDO_SIGMA = 0.005`). Ratios of exactly 0
are common and meaningful for inter-cluster pairs under `heterogeneous`; `scale_y_log10()`
drops them silently, which is what the pseudo-log axis exists to prevent. Do not switch it
back to a plain log scale.

The sliding-window figures (`nrr_by_pair`, `nrr_intra_vs_inter`) are drawn over
`PLOT_TIME_MIN`..`PLOT_TIME_MAX` (90..300) only. The right edge matters: `MAX_TIME = 350`
is a hard stop, and `calc_rr()` needs *both* members of a pair inside the window, so tips
near the end lose the partners that would have completed their pairs and nRR falls for
reasons that have nothing to do with connectivity. The filter is applied to the data
before `geom_smooth()`, not via `xlim`/`coord_cartesian` -- cropping the view would leave
the censored tail still pulling the loess fit down inside the visible range. These bounds
are display-only constants and do not track `MAX_TIME` or `BASELINE_TPOINT`; change them by
hand if either moves.

### Plotting Scripts

- `scripts/<scenario>/*_plot.R` scripts generate summary visualizations after batch runs
- `scripts/two_cluster/two_cluster_combined_plot.R` - stacks all three arms' recovery
  panels into one 6x8 PDF (`recovery_by_phase_combined.pdf`), and likewise the R(t)
  panels (`rt_combined.pdf`). Run after every scenario has been plotted. Scenario display
  names live in `SCENARIO_LABELS` (`none` -> "No Shutdown", `homogeneous` -> "Uniform
  Recovery", `heterogeneous` -> "Staggered Recovery"); the same map is duplicated in
  `two_cluster_plot.R`, so a rename must touch both. `PLOT_TIME_MAX`, `RT_SPAN`, `RT_SEEDING_END`
  (50) and the fixed R(t) axis `RT_LIM` (0..3) are also duplicated between the two
  scripts. R(t) figures start at t = 0 with t < 50 shaded light orange as "Seeding" (not
  "burn-in", which names the t < 100 phase). The seeding spike (R ~ 10) runs off the top
  of the axis via `coord_cartesian()` but stays in the loess fit, which spans the full
  range by choice -- this widens the confidence band across the whole series.
- `scripts/two_cluster/two_cluster_reproduction.R` - per-seed effective reproduction
  number, run by the batch after `_RR.R` (same `--seed` / `--scenario` args). R(t) is the
  mean offspring count of hosts infected in a centred `RT_WINDOW`-day window (7), taken
  from the full `meta.tsv`, so it is independent of subsampling. Hosts with zero
  offspring must be counted -- tallying `inf.by` alone drops them and inflates R(t).
  Appends to `<scenario>_rt.tsv`; the tail near `MAX_TIME` is censored low, so plots crop
  at `PLOT_TIME_MAX`.
- `scripts/two_cluster/two_cluster_network_plot.R` - movement-structure diagram, written
  to `movement_network.pdf` at the repository root rather than under `results/`, since it
  is a property of the model rather than of any sweep. It re-declares `BASE_WEIGHTS`
  because sourcing `two_cluster_sim.R` would run a simulation -- **if that matrix changes
  in the sim, change it here too.**
- `scripts/common/sim_analysis.R` - General purpose analysis utilities
- `scripts/common/plot_pair_time.R` - Temporal analysis of sequence pairs

## External Tools

**phastSim**: Sequence evolution simulator
- Uses `--scale` parameter for molecular clock rate (e.g., 2e-6 for SARS-CoV-2)
- Requires `--ref` FASTA (`data/reference/wuhan.fasta`), `--tree` (Newick), and `--outpath`
- Produces `sars-cov-2_simulation_output.fasta`

**pairsnp**: SNP distance calculator
- `-s`: Output to stdout
- `-d 0`: Filter to identical sequences only
- `-t N`: Number of threads

## Important Notes

- Simulations use `nosoiSim()` with `type="single"` (single-host species)
- Population structure: `"none"` (basic_sim) or `"discrete"` (all structured models)
- The `diff.pTrans`, `timeDep.pTrans`, `hostCount.pTrans` flags control transmission probability complexity
- nosoi validates any `diff.*` callback by grepping its body for the literal text
  `current.in == "<deme>"`, once per deme. Deme lookups in `p_Move_fct` / `p_Trans_fct`
  must therefore stay unrolled as one `if` per deme -- collapsing them into a vector
  index or a `switch()` is a silent validation failure, not a style change.
- Seed management: phastSim typically uses fixed seed (17) while R simulations vary
- Output directories under `results/` are created automatically if they don't exist
- Everything under `results/` is gitignored; `data/` holds only tracked inputs
  (`data/reference/wuhan.fasta`)
