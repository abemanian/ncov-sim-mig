# validate_sampled_tree.R
# Check build_sampled_tree() against nosoi's own tree on small simulations: for a
# random subsample, the directly built tree must match
# keep.tip(getTransmissionTree(sim), ids) in every tip-to-tip distance and every
# root-to-tip depth. Run from the repository root:
#   Rscript scripts/common/validate_sampled_tree.R

suppressPackageStartupMessages({
  library(nosoi, quietly = TRUE)
  library(ape, quietly = TRUE)
})
source("scripts/common/build_sampled_tree.R")

run_small_sim <- function(seed, length_sim, n_init = 1) {
  set.seed(seed)
  nosoiSim(
    type = "single", popStructure = "discrete",
    structure.matrix = matrix(c(0, 1, 1, 0), 2, 2, dimnames = list(c("A", "B"), c("A", "B"))),
    length.sim = length_sim, max.infected = 4000, init.individuals = n_init, init.structure = "A",
    nContact = function(t) rnbinom(1, size = 1, mu = 3), param.nContact = NA,
    pExit = function(t) 0.15, param.pExit = NA,
    pMove = function(t) 0.05, param.pMove = NA,
    pTrans = function(t, t_incub) ifelse(t < t_incub, 0, 0.2),
    param.pTrans = list(t_incub = function(n) rgamma(n, 36, 12)),
    print.progress = FALSE
  )
}

depths <- function(tr) {
  d <- node.depth.edgelength(tr)[seq_len(Ntip(tr))] + (if (is.null(tr$root.edge)) 0 else tr$root.edge)
  setNames(d, tr$tip.label)
}

check <- function(seed, length_sim, n_sample, min_hosts = 500) {
  # Most single-introduction outbreaks die out early; step through seeds until
  # one is big enough to exercise deep chains and many same-day events
  repeat {
    sim <- run_small_sim(seed, length_sim)
    if (getTableHosts(sim)[, .N] >= min_hosts) break
    seed <- seed + 100
  }
  # Deep copy: getTransmissionTree() reorders nosoi's host table by reference,
  # and as.data.frame() alone would still share its column memory
  hosts <- as.data.frame(data.table::copy(getTableHosts(sim)))
  hosts <- hosts[, c("hosts.ID", "inf.by", "inf.time", "out.time")]
  full <- getTransmissionTree(sim)@phylo

  set.seed(seed + 1000)
  ids <- sample(hosts$hosts.ID, min(n_sample, nrow(hosts)))
  ref <- keep.tip(full, ids)
  new <- build_sampled_tree(hosts, ids, sim$total.time)

  dr <- cophenetic(ref)[ids, ids]
  dn <- cophenetic(new)[ids, ids]
  ok_dist <- isTRUE(all.equal(dr, dn, tolerance = 1e-9))
  # Root-to-tip depth (with root edge) must equal each tip's absolute sampling
  # time: out.time, or total.time for hosts still active at the end. Compared to
  # the host table rather than to `ref`, because keep.tip() drops the stem from
  # nosoi's root down to the subsample's common ancestor.
  tip_time <- setNames(hosts$out.time, hosts$hosts.ID)[ids]
  tip_time[is.na(tip_time)] <- sim$total.time
  ok_depth <- isTRUE(all.equal(unname(depths(new)[ids]), unname(tip_time), tolerance = 1e-9))
  ok_shape <- isTRUE(all.equal.phylo(unroot(ref), unroot(new), use.edge.length = TRUE))
  message(sprintf(
    "seed %d: %d hosts, %d sampled, still active %d | distances %s, depths %s, topology %s",
    seed, nrow(hosts), length(ids), sum(is.na(hosts$out.time)),
    ok_dist, ok_depth, ok_shape
  ))
  ok_dist && ok_depth
}

# Several introductions (init.individuals > 1). nosoi has no single tree to
# compare against, so check the two things that define the joined tree: within
# an introduction, distances equal those of that introduction built on its own
# (the single-index case validated above); across introductions, lineages meet
# only at the common root at t0, so d(i, j) = (t_i - t0) + (t_j - t0).
check_multi <- function(seed, length_sim, n_sample, n_init = 5, min_hosts = 500) {
  repeat {
    sim <- run_small_sim(seed, length_sim, n_init)
    if (getTableHosts(sim)[, .N] >= min_hosts) break
    seed <- seed + 100
  }
  hosts <- as.data.frame(data.table::copy(getTableHosts(sim)))[
    , c("hosts.ID", "inf.by", "inf.time", "out.time")]

  # Introduction each host descends from
  intro <- ifelse(hosts$inf.by %in% hosts$hosts.ID, NA, hosts$hosts.ID)
  parent <- match(hosts$inf.by, hosts$hosts.ID)
  while (anyNA(intro)) {
    todo <- which(is.na(intro))
    intro[todo] <- intro[parent[todo]]
  }
  names(intro) <- hosts$hosts.ID

  set.seed(seed + 1000)
  ids <- sample(hosts$hosts.ID, min(n_sample, nrow(hosts)))
  joined <- build_sampled_tree(hosts, ids, sim$total.time)
  dj <- cophenetic(joined)[ids, ids]

  tip_time <- setNames(hosts$out.time, hosts$hosts.ID)
  tip_time[is.na(tip_time)] <- sim$total.time
  t0 <- min(hosts$inf.time[!hosts$inf.by %in% hosts$hosts.ID])

  ok_within <- TRUE
  for (k in unique(intro[ids])) {
    sub_ids <- ids[intro[ids] == k]
    if (length(sub_ids) < 2) next
    sub <- build_sampled_tree(hosts[intro == k, ], sub_ids, sim$total.time)
    ok_within <- ok_within &&
      isTRUE(all.equal(cophenetic(sub)[sub_ids, sub_ids], dj[sub_ids, sub_ids], tolerance = 1e-9))
  }
  cross <- outer(intro[ids], intro[ids], "!=")
  expected <- outer(tip_time[ids] - t0, tip_time[ids] - t0, "+")
  ok_cross <- isTRUE(all.equal(dj[cross], expected[cross], tolerance = 1e-9))
  ok_depth <- isTRUE(all.equal(unname(depths(joined)[ids]), unname(tip_time[ids]), tolerance = 1e-9))
  message(sprintf(
    "seed %d (%d introductions): %d hosts, %d sampled from %d introductions | within %s, across %s, depths %s",
    seed, n_init, nrow(hosts), length(ids), length(unique(intro[ids])),
    ok_within, ok_cross, ok_depth
  ))
  ok_within && ok_cross && ok_depth
}

# Partial subsamples of several sizes, then every host sampled. Runs that end
# at length.sim leave hosts still active, exercising the total.time tip rule.
results <- c(
  check(1, 60, 50), check(2, 60, 300), check(3, 60, 1500),
  check(4, 80, 3000), check(5, 60, 100000),
  check_multi(6, 50, 400), check_multi(7, 50, 2000), check_multi(8, 40, 100000)
)
results <- unlist(results)
if (length(results) == 0 || !all(results)) stop("Validation FAILED")
message("All ", length(results), " checks passed")
