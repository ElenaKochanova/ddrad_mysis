# ============================================================
# run snmf pipeline
#
# Combined pipeline for sNMF analysis:
#   1. Convert a VCF file to GENO format
#   2. Run sNMF for a range of K values
#   3. Read Q-matrix results in R with pophelper
#   4. Plot the results
#   
#
# Requirements:
#   - bash
#   - R with package: pophelper
#   - vcf2geno binary
#   - sNMF binary
#
# Example:
#   bash scripts/run_snmf_pipeline.sh sal_rel.snps.vcf 1 10 1000 2 3
#
# Arguments:
#   $1 = input VCF file
#   $2 = minimum K value for sNMF
#   $3 = maximum K value for sNMF
#   $4 = number of sNMF iterations
#   $5 = K value to plot
#   $6 = number of replicate Q files to plot
#
# Example meaning:
#   sal_rel.snps.vcf 1 10 1000 2 3
#   -> run sNMF for K = 1 to 10, with 1000 iterations,
#      then plot K = 2 using the first 3 replicate files
#
# Notes:
#   - This script assumes the binaries are in ./bin/
#   - It assumes sNMF creates a directory ending in .snmf
#   - Adjust paths if your setup differs
# ============================================================

set -euo pipefail

# ----------------------------
# Input arguments
# ----------------------------
VCF_FILE="${1:-sal_rel.snps.vcf}"
K_MIN="${2:-1}"
K_MAX="${3:-10}"
N_ITER="${4:-1000}"
K_PLOT="${5:-2}"
N_REPS="${6:-3}"

# ----------------------------
# Paths to binaries
# ----------------------------
VCF2GENO="./bin/vcf2geno"
SNMF="./bin/sNMF"

# ----------------------------
# Check required input
# ----------------------------
if [[ ! -f "$VCF_FILE" ]]; then
    echo "ERROR: Input VCF file not found: $VCF_FILE" >&2
    exit 1
fi

if [[ ! -x "$VCF2GENO" ]]; then
    echo "ERROR: vcf2geno binary not found or not executable: $VCF2GENO" >&2
    exit 1
fi

if [[ ! -x "$SNMF" ]]; then
    echo "ERROR: sNMF binary not found or not executable: $SNMF" >&2
    exit 1
fi

# ----------------------------
# Step 1. Convert VCF to GENO
# ----------------------------
echo "Converting VCF to GENO format..."
"$VCF2GENO" "$VCF_FILE"

GENO_FILE="${VCF_FILE%.vcf}.geno"

if [[ ! -f "$GENO_FILE" ]]; then
    echo "ERROR: Expected GENO file was not created: $GENO_FILE" >&2
    exit 1
fi

echo "GENO file created: $GENO_FILE"

# ----------------------------
# Step 2. Run sNMF
# ----------------------------
echo "Running sNMF..."
echo "K range: ${K_MIN}-${K_MAX}"
echo "Iterations: $N_ITER"

"$SNMF" \
    -x "$GENO_FILE" \
    -i "$N_ITER" \
    -K "${K_MIN}-${K_MAX}" \
    -c

echo "sNMF run completed."

# ----------------------------
# Step 3. Define expected result paths
# ----------------------------
# sNMF usually creates an output directory based on the GENO file name.
# Example:
#   sal_rel.snps.geno  -> sal_rel.snps.snmf/
#
# For plotting, we point to:
#   <geno_basename>.snmf/K<K_PLOT>/results
# Adjust this if your sNMF version uses a different folder structure.
# ----------------------------
GENO_BASENAME="$(basename "$GENO_FILE")"
GENO_PREFIX="${GENO_BASENAME%.geno}"
SNMF_DIR="${GENO_PREFIX}.snmf"
RESULTS_DIR="${SNMF_DIR}/K${K_PLOT}/results"
PLOT_FILE="${GENO_PREFIX}_K${K_PLOT}_structure_plot.pdf"

echo "Looking for sNMF results in: $RESULTS_DIR"

if [[ ! -d "$RESULTS_DIR" ]]; then
    echo "ERROR: Results directory not found: $RESULTS_DIR" >&2
    echo "Check the sNMF output structure and adjust RESULTS_DIR in the script if needed." >&2
    exit 1
fi

# Export variables so R can use them
export RESULTS_DIR
export N_REPS
export K_PLOT
export PLOT_FILE


# ============================================================
# step 4. plot the results
#
# Usage:
#   Rscript scripts/plot_snmf.R /path/to/K2/results 3 snmf_K2_plot.pdf
#
# Arguments:
#   1. results directory
#   2. number of replicate files to use
#   3. output PDF filename
# ============================================================

library(pophelper)

args <- commandArgs(trailingOnly = TRUE)

results_dir <- ifelse(length(args) >= 1, args[1],
                      "~/Desktop/Mysis_ddRAD/sNMF/sal_rel_p10_r07.snps.snmf/K2/results")
n_reps <- ifelse(length(args) >= 2, as.integer(args[2]), 3)
outfile <- ifelse(length(args) >= 3, args[3], "snmf_plot.pdf")

sfiles <- list.files(path = results_dir, full.names = TRUE)

if (length(sfiles) == 0) {
  stop("No result files found in: ", results_dir)
}

if (length(sfiles) < n_reps) {
  stop("Requested ", n_reps, " replicates, but only ", length(sfiles), " files found.")
}

slist <- readQ(files = sfiles)
slist_aligned <- alignK(slist[seq_len(n_reps)])

pdf(outfile, width = 10, height = 4)

p1 <- plotQ(
  slist_aligned[seq_len(n_reps)],
  imgoutput = "join",
  showindlab = TRUE,
  returnplot = TRUE,
  exportplot = FALSE,
  basesize = 11,
  sortind = "label",
  sharedindlab = FALSE,
  showyaxis = TRUE,
  showticks = TRUE,
  height = 1.6,
  indlabsize = 2.3,
  indlabheight = 0.8,
  indlabspacer = -1,
  barbordercolour = "white",
  barbordersize = 0,
  clustercol = c("deepskyblue", "yellowgreen")
)

print(p1)
dev.off()
