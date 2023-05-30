# need CI for occ estimates
boot_occ <- function(comm, grp, return_S = FALSE) {
  grp_lvs <- unique(grp)
  row_indices <- NULL
  ct <- NULL
  for (i in seq_along(grp_lvs)) { 
    tmp <- sum(grp == grp_lvs[i])
    print(c(i, tmp))
    ct <- c(ct, tmp)
  }
  for (i in seq_along(grp_lvs)) {
    if (sum(grp == grp_lvs[i]) > 1)
        row_indices <- c(row_indices, sample(which(grp == grp_lvs[i]), replace = TRUE))
    else
        row_indices <- c(row_indices, which(grp == grp_lvs[i]))
  }
  # rand set of comm
  occ <- aggregate(comm[row_indices, ], list(grp), function(x) sum(x > 0) / length(x))
  if (return_S)
    rowSums(occ[ , -1])
  else {
    rownames(occ) <- occ$Group.1
    occ <- occ[ , -1]
    as.matrix(occ)
  }
}

boot_raw_S <- function(comm, grp) {
  grp_lvs <- unique(grp)
  row_indices <- NULL
  for (i in seq_along(grp_lvs)) { 
    row_indices <- c(row_indices, 
                     sample(which(grp == grp_lvs[i]), replace = TRUE))
  }
  # rand set of comm
  S <- rowSums(comm[row_indices, ] > 0)
  Savg <- tapply(S, list(grp), mean)
  Savg
}


boot_sp_occ <- function(comm, grp, return_S = FALSE) {
  grp_lvs <- unique(grp)
  row_indices <- NULL
  for (i in seq_along(grp_lvs)) { 
    row_indices <- c(row_indices, 
                     sample(which(grp == grp_lvs[i]), replace = TRUE))
  }
  # rand set of comm
  occ <- tapply(comm[row_indices], list(grp), function(x) sum(x > 0) / length(x))
  if (return_S)
    sum(occ)
  else {
    t(occ)
  }
}
