# build_sampled_tree.R
# Build the transmission tree of a subsample of hosts directly from a nosoi host
# table, without ever constructing the full tree.
#
# nosoi's getTransmissionTree() builds one tip per infected host and scales
# badly: a 250k-host outbreak takes the better part of an hour and the result is
# then mostly discarded by subsampling. This script follows the same conventions
# but only for the lineages the sample needs:
#
#   * a host's lineage runs from its infection time to out.time (or total.time
#     if it is still active at the end of the run), and ends in a tip;
#   * every onward transmission splits it, and infectees sharing an infector and
#     an infection time share one node (nosoi's `indNodes` grouping);
#   * the root is the index host's first transmission event.
#
# With several initial hosts (init.individuals > 1) nosoi has one tree per
# introduction. Here they are joined at a single root at the earliest index
# infection time (t = 0), with each introduction's lineage hanging off it: all
# introductions carry the same founding genome, so sequences from different
# introductions differ by mutations along both lineages from t = 0.
#
# Only sampled hosts and their ancestors are kept, then single-child nodes are
# collapsed with their branch lengths summed. The result is the tree
# ape::keep.tip(getTransmissionTree(sim), ids) would give, which
# scripts/common/validate_sampled_tree.R checks on a small simulation.

suppressPackageStartupMessages({
  library(ape, quietly = TRUE)
  library(argparse, quietly = TRUE)
})

build_sampled_tree <- function(hosts, ids, total_time) {
  hosts <- as.data.frame(hosts)[, c("hosts.ID", "inf.by", "inf.time", "out.time")]
  hosts$hosts.ID <- as.character(hosts$hosts.ID)
  hosts$inf.by <- as.character(hosts$inf.by)
  hosts$out.time[is.na(hosts$out.time)] <- total_time

  missing <- setdiff(ids, hosts$hosts.ID)
  if (length(missing) > 0) {
    stop(length(missing), " sampled IDs are not in the host table (e.g. ",
         paste(head(missing, 3), collapse = ", "), ")")
  }
  ids <- unique(ids)

  row_of <- setNames(seq_len(nrow(hosts)), hosts$hosts.ID)
  parent_row <- unname(row_of[hosts$inf.by]) # NA for index hosts
  index_rows <- which(is.na(parent_row))
  if (length(index_rows) == 0) {
    stop("No index host (inf.by outside the host table) found")
  }
  t0 <- min(hosts$inf.time[index_rows])

  # Sampled hosts plus all their ancestors, walking inf.by back to the root
  needed <- logical(nrow(hosts))
  frontier <- unname(row_of[ids])
  while (length(frontier) > 0) {
    frontier <- frontier[!needed[frontier]]
    needed[frontier] <- TRUE
    frontier <- unique(parent_row[frontier])
    frontier <- frontier[!is.na(frontier)]
  }
  need_rows <- which(needed)

  # Tips first (1..n), as ape requires, then one node per transmission event
  n_tip <- length(ids)
  tip_of_row <- integer(nrow(hosts))
  tip_of_row[row_of[ids]] <- seq_len(n_tip)

  # Event = (infector, infection time) among the needed non-index hosts
  child_rows <- setdiff(need_rows, index_rows)
  ev_key <- paste(parent_row[child_rows], hosts$inf.time[child_rows], sep = "@")
  ev_levels <- unique(ev_key)
  ev_id <- match(ev_key, ev_levels) # per child row
  ev_host <- parent_row[child_rows][match(ev_levels, ev_key)]
  ev_time <- hosts$inf.time[child_rows][match(ev_levels, ev_key)]
  n_ev <- length(ev_levels)
  ev_node <- n_tip + seq_len(n_ev) # provisional node numbers

  # One extra node: the common root joining every introduction, at t0
  root <- n_tip + n_ev + 1L
  node_time <- c(hosts$out.time[row_of[ids]], ev_time, t0)

  # Each needed host's chain: its own events in time order, then its tip if
  # sampled. Consecutive elements are joined by an edge.
  ord <- order(ev_host, ev_time)
  chain_host <- ev_host[ord]
  chain_node <- ev_node[ord]
  chain_time <- ev_time[ord]
  sampled_rows <- unname(row_of[ids])
  chain_host <- c(chain_host, sampled_rows)
  chain_node <- c(chain_node, tip_of_row[sampled_rows])
  chain_time <- c(chain_time, hosts$out.time[sampled_rows])
  # A host can transmit and exit in the same day, so on a time tie the tip must
  # sort after the event: a tip with a child edge is a malformed tree that
  # crashes ape's C routines rather than raising an error
  ord <- order(chain_host, chain_time, chain_node <= n_tip)
  chain_host <- chain_host[ord]
  chain_node <- chain_node[ord]

  same_host <- chain_host[-1] == chain_host[-length(chain_host)]
  e_chain <- cbind(chain_node[-length(chain_node)][same_host], chain_node[-1][same_host])

  # Each needed child's lineage enters at its infection event, and its first
  # chain element hangs off that event node
  first_idx <- !duplicated(chain_host)
  first_node_of_row <- integer(nrow(hosts))
  first_node_of_row[chain_host[first_idx]] <- chain_node[first_idx]
  e_entry <- cbind(ev_node[ev_id], first_node_of_row[child_rows])

  # Each introduction that has sampled descendants hangs off the common root.
  # With one such introduction the root is a single-child node, which
  # collapse.singles() below folds away.
  needed_index <- index_rows[needed[index_rows]]
  e_root <- cbind(root, first_node_of_row[needed_index])

  edge <- rbind(e_chain, e_entry, e_root)

  # Renumber internal nodes so the root is n_tip + 1, as ape expects
  internal <- c(root, ev_node)
  remap <- c(seq_len(n_tip), integer(n_ev + 1L))
  remap[internal] <- n_tip + seq_along(internal)
  edge_len <- node_time[edge[, 2]] - node_time[edge[, 1]]

  tree <- structure(list(
    edge = matrix(remap[edge], ncol = 2),
    edge.length = edge_len,
    tip.label = ids,
    Nnode = n_ev + 1L,
    root.edge = 0
  ), class = "phylo")
  if (any(tree$edge[, 1] <= n_tip) || anyDuplicated(tree$edge[, 2]) ||
      any(edge_len < 0)) {
    stop("build_sampled_tree: malformed tree (tip as parent, node with two ",
         "parents, or negative branch length)")
  }
  tree <- reorder(tree, "cladewise")
  tree <- collapse.singles(tree, root.edge = TRUE)
  # collapse.singles() discards an existing root edge when it collapses a
  # single-child root, so recompute it: the stem runs from t0 to the sample's
  # common ancestor, making every root-to-tip depth equal to that tip's
  # absolute sampling time (and zero when several introductions are joined)
  tree$root.edge <- (node_time[1] - t0) - node.depth.edgelength(tree)[1]
  tree
}

if (sys.nframe() == 0) {
  parser <- ArgumentParser()
  parser$add_argument("--meta", required = TRUE, help = "nosoi host table (meta.tsv)")
  parser$add_argument("--ids", required = TRUE, help = "Host IDs to keep, one per line")
  parser$add_argument("--total-time", type = "double", required = TRUE,
                      help = "Simulation end time; out.time of still-active hosts")
  parser$add_argument("--output", required = TRUE, help = "Output Newick tree")
  args <- parser$parse_args()

  meta <- read.delim(args$meta, stringsAsFactors = FALSE)
  ids <- readLines(args$ids)
  tree <- build_sampled_tree(meta, ids, args$total_time)
  write.tree(tree, args$output)
  message("build_sampled_tree: ", Ntip(tree), " tips, ", tree$Nnode,
          " internal nodes, from ", nrow(meta), " hosts")
}
