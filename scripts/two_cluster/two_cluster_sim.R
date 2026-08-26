options(tidyverse.quiet = TRUE)
qlibrary <- function(pkg, ...) {
  suppressPackageStartupMessages(library(pkg,
  quietly = TRUE,
  character.only = TRUE, ...))
}

qlibrary("nosoi")
qlibrary("tidyverse")
qlibrary("viridis")
qlibrary("ape")
qlibrary("tidytree")
qlibrary("treeio")
qlibrary("argparse")

collect_args <- function() {
  parser <- ArgumentParser()
  parser$add_argument("--seed",
    type = "integer",
    default = 66,
    help = "Simulation Seed")
  parser$add_argument("--draw",
    action = "store_true",
    default = FALSE,
    help = "Draw plots")
  parser$add_argument("--scenario",
    type = "character",
    default = "homogeneous",
    help = "Scenario: none, homogeneous, heterogeneous")
  parser$parse_args()
}

args <- collect_args()
sim_seed <- args$seed
draw_plots <- args$draw
scenario <- args$scenario

VALID_SCENARIOS <- c("none", "homogeneous", "heterogeneous")
if (!scenario %in% VALID_SCENARIOS) {
  stop("--scenario must be one of: ", paste(VALID_SCENARIOS, collapse = ", "))
}

# Display name for plot titles; `scenario` itself stays as-is for filenames
scenario_label <- function(s) if (s == "none") "No Shutdown" else s

# ---------------------------------------------------------------------------
# Population structure
# ---------------------------------------------------------------------------
# Two fully-connected triangles joined by a single direct edge, A <-> X. There
# is no intermediate bridge deme: A and X are each other's only route between
# clusters, and each is also an ordinary member of its own triangle.
#
#    B----C          Y----Z
#     \  /            \  /
#      A ---- 0.3 ---- X

DEMES <- c("A", "B", "C", "X", "Y", "Z")
CLUSTER_1 <- c("A", "B", "C")
CLUSTER_2 <- c("X", "Y", "Z")
GATEWAYS <- c("A", "X")

CLUSTER_OF <- c(A = "C1", B = "C1", C = "C1", X = "C2", Y = "C2", Z = "C2")

# Baseline probability that a host moves at all, per timestep
PROB_MOVE <- 0.10

# Destination weights conditional on moving (rows sum to 1)
BASE_WEIGHTS <- matrix(
  c(
    0.00, 0.35, 0.35, 0.30, 0.00, 0.00, # A -> B .35  C .35  X .30
    0.50, 0.00, 0.50, 0.00, 0.00, 0.00, # B -> A .50  C .50
    0.50, 0.50, 0.00, 0.00, 0.00, 0.00, # C -> A .50  B .50
    0.30, 0.00, 0.00, 0.00, 0.35, 0.35, # X -> A .30  Y .35  Z .35
    0.00, 0.00, 0.00, 0.50, 0.00, 0.50, # Y -> X .50  Z .50
    0.00, 0.00, 0.00, 0.50, 0.50, 0.00  # Z -> X .50  Y .50
  ),
  nrow = 6, byrow = TRUE,
  dimnames = list(DEMES, DEMES)
)

# Per-edge movement RATES: rate_ij = P(move) * P(destination j | move).
# Phases are defined directly on these rates, which is what makes an
# edge-specific policy expressible: scaling one edge changes both that deme's
# total mobility and its destination mix, in the right proportions.
BASE_RATES <- BASE_WEIGHTS * PROB_MOVE

# The A <-> X edge, the only inter-cluster connection
BRIDGE_EDGES <- rbind(c("A", "X"), c("X", "A"))

# ---------------------------------------------------------------------------
# Intervention schedule
# ---------------------------------------------------------------------------
# Five phases. Multipliers are fractions of the baseline per-edge rate.
# Suppressed travel is LOST, not redirected: restricting the A-X edge lowers
# A's total mobility and leaves the absolute A->B and A->C rates untouched, so
# any change on within-cluster edges is not an artefact of the bridge policy.
PHASE_STARTS <- c(0, 100, 150, 200, 250)
PHASE_NAMES <- c("burn-in", "shutdown", "initial", "moderate", "full")

# Applied to every edge
PHASE_ALL <- c(1.00, 0.10, 0.33, 0.66, 1.00)
# Applied to the A-X edge under the heterogeneous scenario: restoration lags
# a full phase behind, and never gets past halfway.
PHASE_BRIDGE_HETERO <- c(1.00, 0.10, 0.10, 0.33, 0.50)

phase_index <- function(prestime) findInterval(prestime, PHASE_STARTS)

# Per-edge rate matrix in force at a given time
phase_rates <- function(prestime) {
  if (scenario == "none") {
    return(BASE_RATES)
  }
  k <- phase_index(prestime)
  rates <- BASE_RATES * PHASE_ALL[k]
  if (scenario == "heterogeneous") {
    rates[BRIDGE_EDGES] <- BASE_RATES[BRIDGE_EDGES] * PHASE_BRIDGE_HETERO[k]
  }
  rates
}

# Probability of moving at all, by deme, at a given time
phase_pmove <- function(prestime) rowSums(phase_rates(prestime))

# nosoi tests structure.matrix row sums with `!= 1`, so they must be exactly 1
# in floating point. Dividing by the row total is not enough: a row such as
# 0.30 + 0.35 + 0.35 lands one ULP short and stays short after normalising, so
# the residual slack is folded back into each row's largest entry.
exact_rows <- function(m) {
  for (i in seq_len(nrow(m))) {
    row <- m[i, ]
    for (k in 1:10) {
      slack <- 1 - sum(row)
      if (slack == 0) break
      j <- which.max(row)
      row[j] <- row[j] + slack
    }
    m[i, ] <- row
  }
  m
}

# Destination weights, by deme, at a given time
phase_matrix <- function(prestime) {
  rates <- phase_rates(prestime)
  total <- rowSums(rates)
  # A deme with no outbound movement still needs a valid row; nosoi only ever
  # samples it for hosts that already drew a move, which cannot happen at rate 0
  total[total == 0] <- 1
  exact_rows(rates / total)
}

# ---------------------------------------------------------------------------
# Time-varying movement: patch nosoi's discrete move function
# ---------------------------------------------------------------------------
# nosoi fixes structure.matrix for the whole run, so destination probabilities
# are normally constant. moveFunction.discrete does read the matrix afresh on
# every call and has pres.time in scope, so injecting a single line at the top
# of its body makes the matrix time-varying while leaving the rest of nosoi's
# logic untouched. The original body is spliced in programmatically rather than
# retyped, so this cannot drift from the installed source.
patch_move_function <- function() {
  orig <- nosoi:::moveFunction.discrete
  required <- c("res", "pres.time", "moving.full", "structure.matrix")
  if (!all(required %in% names(formals(orig)))) {
    stop(
      "nosoi:::moveFunction.discrete has an unexpected signature; ",
      "the time-varying movement patch needs updating for nosoi ",
      as.character(utils::packageVersion("nosoi"))
    )
  }
  patched <- orig
  body(patched) <- as.call(c(
    as.name("{"),
    quote(structure.matrix <- phase_matrix(pres.time)),
    as.list(body(orig))[-1]
  ))
  environment(patched) <- new.env(parent = environment(orig))
  assign("phase_matrix", phase_matrix, envir = environment(patched))
  utils::assignInNamespace("moveFunction.discrete", patched, ns = "nosoi")
  invisible(TRUE)
}

patch_move_function()

# ---------------------------------------------------------------------------
# Epidemiological constants
# ---------------------------------------------------------------------------
INCUB_MEAN <- 3
INCUB_SD <- 0.5
PROB_EXIT <- 1.0 / 7.0
CONTACT_MEAN <- 4
CONTACT_K <- 0.2

INF_ALPHA <- 1
INF_BETA <- 4

MAX_CAPACITY <- 300
CAPACITY_MIDPOINT <- MAX_CAPACITY / 2
CAPACITY_SCALE <- MAX_CAPACITY / 10

calc_capacity_factor <- function(n) {
  1 / (1 + exp((n - CAPACITY_MIDPOINT) / CAPACITY_SCALE))
}

# ---------------------------------------------------------------------------
# nosoi callbacks
# ---------------------------------------------------------------------------
# NOTE: nosoi validates any diff.* function by grepping its body for the literal
# text `current.in == "<deme>"` for every deme, hence the unrolled lookups.

p_Exit_fct <- function(t) {
  PROB_EXIT
}

p_Move_fct <- function(t, prestime, current.in) {
  pmove <- phase_pmove(prestime)

  if (current.in == "A") p <- pmove[["A"]]
  if (current.in == "B") p <- pmove[["B"]]
  if (current.in == "C") p <- pmove[["C"]]
  if (current.in == "X") p <- pmove[["X"]]
  if (current.in == "Y") p <- pmove[["Y"]]
  if (current.in == "Z") p <- pmove[["Z"]]

  return(p)
}

n_Contact_fct <- function(t) rnbinom(1, size = CONTACT_K, mu = CONTACT_MEAN)

p_Trans_fct <- function(t, prestime, current.in, host.count, t_incub) {
  # Movement is the only thing the intervention touches; transmission is
  # identical in every deme and every phase.
  if (current.in == "A") place_factor <- 1
  if (current.in == "B") place_factor <- 1
  if (current.in == "C") place_factor <- 1
  if (current.in == "X") place_factor <- 1
  if (current.in == "Y") place_factor <- 1
  if (current.in == "Z") place_factor <- 1

  base_p <- rbeta(1, INF_ALPHA, INF_BETA)
  capacity_factor <- calc_capacity_factor(host.count)
  p <- base_p * capacity_factor * place_factor
  ifelse(t < t_incub, 0, p)
}

d_Incub <- function(n, mean = INCUB_MEAN, sd = INCUB_SD) {
  shape <- (mean / sd)^2
  rate <- shape / mean
  rgamma(n, shape = shape, rate = rate)
}

param_pTrans <- list(t_incub = d_Incub)

# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------
set.seed(sim_seed)
MAX_OUTBREAK <- 120000
MAX_TIME <- 350
MAX_TRIES <- 50
MIN_INFECTED <- 200
# Acceptance is judged on the burn-in period ONLY. Every scenario is identical
# before T_INTERVENTION, so all three accept or reject exactly the same
# attempts and end up genuinely paired by seed. Testing anything measured after
# the intervention starts (e.g. final cluster-2 size) would let the more
# restrictive scenarios reject more often and silently accept a different
# attempt, with a different burn-in, breaking the cross-scenario comparison.
T_INTERVENTION <- PHASE_STARTS[2]
MIN_INFECTED_BURN_IN <- 100

INIT_DEME <- "A"

simulate_until_size <- function(tries, min_inf, min_burn_in) {
  for (k in 1:tries) {
    sim <- nosoiSim(
      type = "single",
      popStructure = "discrete",
      structure.matrix = phase_matrix(0),
      length.sim = MAX_TIME,
      max.infected = MAX_OUTBREAK,
      init.individuals = 1,
      init.structure = INIT_DEME,

      nContact = n_Contact_fct,
      param.nContact = NA,
      timeDep.nContact = FALSE,

      pExit = p_Exit_fct,
      param.pExit = NA,
      timeDep.pExit = FALSE,

      pTrans = p_Trans_fct,
      param.pTrans = param_pTrans,
      timeDep.pTrans = TRUE,
      diff.pTrans = TRUE,
      hostCount.pTrans = TRUE,

      pMove = p_Move_fct,
      param.pMove = NA,
      timeDep.pMove = TRUE,
      diff.pMove = TRUE,

      prefix.host = "H",
      print.progress = TRUE
    )
    hosts <- getTableHosts(sim)
    burn_in_size <- sum(hosts$inf.time < T_INTERVENTION)
    if (nrow(hosts) >= min_inf && burn_in_size >= min_burn_in) {
      return(sim)
    }
  }
  stop("Failed to reach target size after ", tries, " attempts.")
}

two_cluster_sim <- simulate_until_size(
  tries = MAX_TRIES,
  min_inf = MIN_INFECTED,
  min_burn_in = MIN_INFECTED_BURN_IN
)

# getTableHosts() returns nosoi's internal data.table by reference and
# getTransmissionTree() reorders it in place, so build the tree first and
# detach with as_tibble() before deriving any column.
tree_sim <- getTransmissionTree(two_cluster_sim)

data_sim <- getTableHosts(two_cluster_sim) %>%
  as_tibble() %>%
  mutate(cluster = unname(CLUSTER_OF[as.character(inf.in)]))

stopifnot(identical(
  data_sim$cluster,
  unname(CLUSTER_OF[as.character(data_sim$inf.in)])
))

# Realised schedule, for the plotting stage
schedule_log <- tibble(
  phase = PHASE_NAMES,
  start = PHASE_STARTS,
  end = c(PHASE_STARTS[-1], MAX_TIME),
  within_mult = if (scenario == "none") 1 else PHASE_ALL,
  bridge_mult = if (scenario == "none") {
    1
  } else if (scenario == "heterogeneous") {
    PHASE_BRIDGE_HETERO
  } else {
    PHASE_ALL
  }
) %>%
  mutate(seed = sim_seed, scenario = scenario)

if (!dir.exists("results/two_cluster")) {
  dir.create("results/two_cluster", recursive = TRUE)
}

write.tree(tree_sim@phylo, "results/two_cluster/tree.nwk")
write_tsv(data_sim, "results/two_cluster/meta.tsv")
write_tsv(schedule_log, "results/two_cluster/schedule.tsv")

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
if (draw_plots) {
  marks <- tibble(t = PHASE_STARTS[-1])

  data_sim %>%
    count(inf.time, cluster, name = "n") %>%
    ggplot(aes(x = inf.time, y = n, color = cluster)) +
    geom_vline(
      data = marks, aes(xintercept = t),
      inherit.aes = FALSE, linetype = "dashed", color = "grey40"
    ) +
    geom_line(alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.85) +
    labs(
      x = "Simulation Time", y = "New Cases", color = "Cluster",
      title = paste0("Incidence by cluster (scenario: ", scenario_label(scenario), ")")
    ) +
    theme_bw()
  ggsave(
    filename = paste0("results/two_cluster/", scenario, "_incidence.jpg"),
    width = 10, height = 6, dpi = 192
  )

  data_sim %>%
    count(inf.time, inf.in, name = "n") %>%
    ggplot(aes(x = inf.time, y = n)) +
    geom_vline(
      data = marks, aes(xintercept = t),
      linetype = "dashed", color = "grey40"
    ) +
    geom_line(alpha = 0.8, color = "firebrick") +
    facet_wrap(vars(inf.in), ncol = 2, scales = "free_y") +
    labs(
      x = "Simulation Time", y = "New Cases",
      title = paste0("Incidence by deme (scenario: ", scenario_label(scenario), ")")
    ) +
    theme_bw()
  ggsave(
    filename = paste0("results/two_cluster/", scenario, "_incidence_deme.jpg"),
    width = 10, height = 8, dpi = 192
  )

  # The policy itself: per-edge movement rate over time
  schedule_log %>%
    select(phase, start, end, within_mult, bridge_mult) %>%
    pivot_longer(c(within_mult, bridge_mult), names_to = "edge", values_to = "mult") %>%
    mutate(edge = if_else(edge == "within_mult", "Within-cluster", "A-X bridge")) %>%
    ggplot(aes(x = start, xend = end, y = mult, yend = mult, color = edge)) +
    geom_segment(linewidth = 2) +
    geom_vline(
      data = marks, aes(xintercept = t),
      inherit.aes = FALSE, linetype = "dashed", color = "grey40"
    ) +
    scale_color_manual(values = c("Within-cluster" = "#38c", "A-X bridge" = "#e94")) +
    ylim(0, 1.1) +
    labs(
      x = "Simulation Time", y = "Fraction of baseline movement rate",
      color = NULL,
      title = paste0("Movement schedule (scenario: ", scenario_label(scenario), ")")
    ) +
    theme_bw()
  ggsave(
    filename = paste0("results/two_cluster/", scenario, "_schedule.jpg"),
    width = 10, height = 5, dpi = 192
  )
}
