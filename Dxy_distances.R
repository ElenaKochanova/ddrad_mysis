# Calculate pairwise individual distances, Dxy, and Dst from a VCF
# for ddRAD data.
#
# Input files:
#   - mysis4sp_r08.vcf.gz
#   - popmap_4sp.txt   (2 columns: sample, population)
#
# Output files:
#   - Dxy_ddRAD.txt
#   - Dst_ddRAD.txt
#   - Dxy_pct_ddRAD.txt
#   - Dst_pct_ddRAD.txt

# Install if needed:
# install.packages("vcfR")

library(vcfR)

# -----------------------------
# 1. Read VCF and extract GTs
# -----------------------------
vcf <- read.vcfR("mysis4sp_r08.vcf.gz")

# Extract genotype strings such as 0/0, 0/1, 1/1, ./.
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

# Convert diploid GT strings to alternate allele counts:
# 0/0 -> 0
# 0/1 or 1/0 -> 1
# 1/1 -> 2
# missing -> NA
gt_to_altcount <- function(x) {
  x <- sub(":.*$", "", x)   # remove anything after ":" if present
  x[x %in% c(".", "./.", ".|.")] <- NA
  
  out <- rep(NA_real_, length(x))
  out[x %in% c("0/0", "0|0")] <- 0
  out[x %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
  out[x %in% c("1/1", "1|1")] <- 2
  
  return(out)
}

geno <- apply(gt, 2, gt_to_altcount)
rownames(geno) <- rownames(gt)

# Samples from VCF
samples <- colnames(geno)

# -----------------------------
# 2. Read and validate popmap
# -----------------------------
popmap <- read.table(
  "popmap_4sp.txt",
  header = FALSE,
  stringsAsFactors = FALSE,
  fill = TRUE
)

colnames(popmap) <- c("sample", "population")

# Basic checks
stopifnot(ncol(popmap) == 2)
stopifnot(all(samples %in% popmap$sample))

# Reorder popmap to match VCF sample order
popmap <- popmap[match(samples, popmap$sample), ]

# Confirm alignment
stopifnot(identical(samples, popmap$sample))

cat("Population counts:\n")
print(table(popmap$population))

# -----------------------------
# 3. Pairwise individual distances
# -----------------------------
# Distance per SNP:
# |gi - gj| / 2 gives:
#   0   for same genotype dosage
#   0.5 for one allele difference
#   1   for opposite homozygotes
#
# Final distance = mean across shared non-missing SNPs.

n_ind <- ncol(geno)

ind_dist <- matrix(
  NA_real_, n_ind, n_ind,
  dimnames = list(samples, samples)
)

diag(ind_dist) <- 0

for (i in 1:(n_ind - 1)) {
  for (j in (i + 1):n_ind) {
    gi <- geno[, i]
    gj <- geno[, j]
    
    ok <- !is.na(gi) & !is.na(gj)
    if (!any(ok)) next
    
    d <- abs(gi[ok] - gj[ok]) / 2
    ind_dist[i, j] <- mean(d)
    ind_dist[j, i] <- ind_dist[i, j]
  }
}

cat("\nSummary of pairwise individual distances:\n")
print(summary(ind_dist[upper.tri(ind_dist)]))

# Optional output
write.table(
  ind_dist,
  file = "individual_distance_matrix.txt",
  sep = "\t", quote = FALSE, col.names = NA
)

# -----------------------------
# 4. Allele frequencies by population
# -----------------------------
pop_levels <- unique(popmap$population)

idx_by_pop <- lapply(pop_levels, function(p) {
  which(popmap$population == p)
})
names(idx_by_pop) <- pop_levels

# For each SNP and population:
# p = alternate allele frequency
p_mat <- sapply(idx_by_pop, function(idx) {
  sub <- geno[, idx, drop = FALSE]
  
  # Sum alternate alleles across individuals
  alt <- rowSums(sub, na.rm = TRUE)
  
  # Number of non-missing diploid genotypes
  n <- rowSums(!is.na(sub))
  
  # Alternate allele frequency
  p <- alt / (2 * n)
  p[n == 0] <- NA
  
  p
})

# Ensure matrix format even if only two populations
p_mat <- as.matrix(p_mat)
colnames(p_mat) <- pop_levels

# Within-population nucleotide diversity (pi)
pi_pop <- apply(p_mat, 2, function(p) {
  mean(2 * p * (1 - p), na.rm = TRUE)
})

cat("\nWithin-population pi values:\n")
print(pi_pop)

# -----------------------------
# 5. Dxy and Dst among populations
# -----------------------------
n_pop <- length(pop_levels)

Dxy <- matrix(
  NA_real_, n_pop, n_pop,
  dimnames = list(pop_levels, pop_levels)
)

Dst <- matrix(
  NA_real_, n_pop, n_pop,
  dimnames = list(pop_levels, pop_levels)
)

diag(Dxy) <- NA
diag(Dst) <- 0

for (i in 1:(n_pop - 1)) {
  for (j in (i + 1):n_pop) {
    p1 <- p_mat[, i]
    p2 <- p_mat[, j]
    
    ok <- !is.na(p1) & !is.na(p2)
    if (!any(ok)) next
    
    # Dxy per site:
    # p1(1-p2) + p2(1-p1)
    dxy_sites <- p1[ok] * (1 - p2[ok]) + p2[ok] * (1 - p1[ok])
    Dxy_ij <- mean(dxy_sites)
    
    # Within-population diversity over the same sites
    pi1 <- mean(2 * p1[ok] * (1 - p1[ok]))
    pi2 <- mean(2 * p2[ok] * (1 - p2[ok]))
    
    # Net divergence
    Dst_ij <- Dxy_ij - (pi1 + pi2) / 2
    
    Dxy[i, j] <- Dxy_ij
    Dxy[j, i] <- Dxy_ij
    
    Dst[i, j] <- Dst_ij
    Dst[j, i] <- Dst_ij
  }
}

cat("\nDxy matrix:\n")
print(Dxy)

cat("\nDst matrix:\n")
print(Dst)

# -----------------------------
# 6. Optional percent versions
# -----------------------------
Dxy_pct <- 100 * Dxy
Dst_pct <- 100 * Dst

cat("\nDxy (%):\n")
print(Dxy_pct)

cat("\nDst (%):\n")
print(Dst_pct)

# -----------------------------
# 7. Write outputs
# -----------------------------
write.table(
  Dxy,
  file = "Dxy_ddRAD.txt",
  sep = "\t", quote = FALSE, col.names = NA
)

write.table(
  Dst,
  file = "Dst_ddRAD.txt",
  sep = "\t", quote = FALSE, col.names = NA
)

write.table(
  Dxy_pct,
  file = "Dxy_pct_ddRAD.txt",
  sep = "\t", quote = FALSE, col.names = NA
)

write.table(
  Dst_pct,
  file = "Dst_pct_ddRAD.txt",
  sep = "\t", quote = FALSE, col.names = NA
)
