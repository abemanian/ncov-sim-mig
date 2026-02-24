# ncov-sim-mig

## Methods

### Simulation Framework

We validated the normalized relative risk (nRR) approach using individual-based epidemiological simulations implemented with the `nosoi` package (v1.1.3) in R. Simulations modeled pathogen transmission across a structured metapopulation of three discrete demes (A, O, B), where deme O served as a central hub. To validate nRR's ability to detect changes in connectivity between sub-populations, we compared a baseline scenario against an isolation scenario in which inter-deme movement was experimentally suppressed.

### Population Structure and Movement

Hosts occupied one of three demes and could move between demes at each time step with a baseline probability of `p_move = 0.02`. Conditional on movement, a host's destination was drawn according to the following asymmetric movement matrix:

|        | → A | → O | → B |
|--------|-----|-----|-----|
| **A:** | —   | 0.8 | 0.2 |
| **O:** | 0.5 | —   | 0.5 |
| **B:** | 0.2 | 0.8 | —   |

Deme O therefore acted as a central transit hub with equal connectivity to both A and B, while A and B were predominantly linked through O.

### Epidemiological Parameters

Transmission dynamics included density-dependent suppression via a logistic carrying capacity function. The per-contact transmission probability at each time step was drawn as:

```
p_trans = Beta(α=4, β=8) × capacity_factor(n) × place_scale
```

where `capacity_factor(n) = 1 / (1 + exp((n − 250) / 50))` is a logistic decline in transmissibility as the local active host count `n` approaches the carrying capacity midpoint of 250 (out of a maximum of `MAX_CAPACITY = 500`). Transmission was set to zero during the incubation period. The full parameter set was:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `BASE_INF_ALPHA` | 4 | Beta distribution shape (α) for transmission probability |
| `BASE_INF_BETA` | 8 | Beta distribution shape (β) for transmission probability |
| `MAX_CAPACITY` | 500 | Logistic carrying capacity |
| `CONTACT_MEAN` | 1.5 | Mean contacts per time step (Poisson) |
| `PROB_EXIT` | 1/7 ≈ 0.143 | Per-step probability of recovery/exit (mean infectious period = 7 steps) |
| `INCUB_MEAN` | 3 | Mean incubation period (time steps) |
| `INCUB_SD` | 0.5 | Standard deviation of incubation period |
| `PROB_MOVE` | 0.02 | Per-step movement probability |

The incubation period was drawn from a gamma distribution parameterized by method of moments: shape = `(mean/sd)² = 36`, rate = `shape/mean = 12`.

### Isolation Scenario

To test nRR's sensitivity to abrupt changes in connectivity, we implemented an isolation scenario in which deme A's outward movement probability was set to zero (`ISOLATED_MOVEMENT = 0`) after simulation time step 150 (`ISOLATED_TIME = 150`), while demes O and B retained their baseline movement probability. This created a step-change in the epidemiological linkage of deme A to the broader metapopulation mid-outbreak, providing a ground-truth signal against which nRR could be evaluated. Demes O and B remained fully connected throughout.

### Simulation Execution and Acceptance Criteria

Each simulation was initialized with a single infected host in deme O and run for up to `MAX_TIME = 1000` time steps or `MAX_OUTBREAK = 30,000` total infections. Simulations yielding fewer than `MIN_INFECTED = 50` total hosts were rejected and re-run (up to 50 attempts). Each scenario was replicated across 50 independent random seeds (seeds 1–50) to quantify variability in nRR estimates.

### Sequence Evolution and SNP Analysis

For each simulation replicate, the resulting transmission tree (Newick format) was passed to **phastSim** to simulate nucleotide sequence evolution under a molecular clock, using the SARS-CoV-2 Wuhan-Hu-1 reference genome (`wuhan.fasta`) as a starting sequence. A fixed substitution rate of `2 × 10⁻⁶` substitutions per site per day (calibrated for SARS-CoV-2) was applied via the `--scale 0.000002` flag, and the phastSim seed was held constant at 17 across all replicates to isolate stochasticity to the epidemiological component.

Pairwise genomic distances between all sampled sequences were computed using **pairsnp**. Only sequence pairs with zero SNP distance (i.e., genetically identical pairs, `-d 0`) were retained, as these are the most informative signal for recent shared transmission history.

### Relative Risk and nRR Calculation

From the identical-sequence pairs, we computed a geographic relative risk (RR) for each ordered deme pair (X, Y):

```
RR(X, Y) = (n_XY × N) / (n_X × n_Y)
```

where `n_XY` is the number of identical pairs with one sequence from deme X and one from deme Y, `n_X` and `n_Y` are the total identical pairs involving deme X and Y respectively, and `N` is the total count across all deme pairs. RR > 1 indicates more genetic similarity between demes X and Y than expected under random mixing.

To control for variation in within-deme transmission intensity—which can inflate or deflate raw RR values independently of between-deme connectivity—we computed the normalized RR (nRR):

```
nRR(X, Y) = RR(X, Y) / mean(RR(X, X), RR(Y, Y))
```

nRR is thus normalized by the arithmetic mean of the two relevant within-deme RR values. A sliding time window (width ±20 time steps, advancing every 5 time steps) was applied to compute nRR as a function of simulation time, enabling detection of the step-change in deme A connectivity at time step 150.

### Software

All simulations were implemented in **R** (≥4.0) using the following packages: `nosoi` (individual-based transmission simulation), `ape` and `treeio`/`ggtree` (phylogenetic tree handling and visualization), `tidyverse` (data manipulation), `igraph`/`ggraph` (transmission network plotting), and `argparse` (command-line argument parsing). Sequence evolution was simulated with **phastSim** and pairwise SNP distances were calculated with **pairsnp**. Analysis scripts were executed in a reproducible `conda` environment (`ncov-sim`) using Python 3.12.
