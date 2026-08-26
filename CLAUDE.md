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

All simulations use the `nosoi` package to model individual-based transmission dynamics with three spatial demes (A, O, B) connected by a movement matrix.

1. **scripts/basic/basic_sim.R** - Single-population baseline simulation with no spatial structure
2. **scripts/shutdown_threebody/shutdown_threebody_sim.R** - Models intervention scenarios with time-varying contact rates and movement restrictions
3. **scripts/carrying_capacity/carrying_capacity_sim.R** - Includes density-dependent transmission via logistic carrying capacity, plus optional seasonal forcing, isolation, and immunization scenarios

### Common Simulation Parameters

All simulations define:
- `p_Exit_fct`: Probability of host exit (death/recovery)
- `n_Contact_fct`: Number of contacts per time step (Poisson distributed)
- `p_Trans_fct`: Transmission probability (often beta-distributed, may depend on incubation period and time)
- `d_Incub`: Incubation period distribution (gamma)
- `p_Move_fct`: Movement probability between demes (three-deme models only)

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
3. **Pairwise SNP analysis** via pairsnp (identify identical/similar sequences)
4. **Relative risk calculation** using `scripts/common/calc_rr.R` (geographic association analysis)

### Batch Execution

Shell scripts run parameter sweeps across multiple seeds:
- `scripts/shutdown_threebody/shutdown_threebody_batch.sh` - Shutdown scenario iterations (50 seeds)
- `scripts/carrying_capacity/carrying_capacity_batch.sh` - Carrying capacity scenario iterations

Pattern: Each script loops over seeds, runs simulation → phastSim → pairsnp → RR calculation → cleanup

### Command Line Arguments

Most R scripts accept:
- `--seed INTEGER`: Random seed (default: 66)
- `--draw`: Generate diagnostic plots (default: FALSE)

Example:
```bash
Rscript scripts/shutdown_threebody/shutdown_threebody_sim.R --seed 123 --draw
```

## Key Functions and Utilities

### scripts/common/calc_rr.R

Calculates geographic relative risk (RR) from pairwise identical sequences:
- `calc_rr()`: Computes RR matrix for sequence pairs between demes
- `normalized_rr()`: Normalizes RR by diagonal entries (nRR metric)
- Supports time window filtering for temporal analysis

RR interpretation: RR(X,Y) > 1 indicates more genetic similarity between demes X and Y than expected by chance

### Plotting Scripts

- `scripts/<scenario>/*_plot.R` scripts generate summary visualizations after batch runs
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
- Population structure: `"none"` (basic_sim) or `"discrete"` (three-deme models)
- The `diff.pTrans`, `timeDep.pTrans`, `hostCount.pTrans` flags control transmission probability complexity
- Seed management: phastSim typically uses fixed seed (17) while R simulations vary
- Output directories under `results/` are created automatically if they don't exist
- Everything under `results/` is gitignored; `data/` holds only tracked inputs
  (`data/reference/wuhan.fasta`)
