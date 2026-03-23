
###############################################################################
# ddRAD de novo assembly pipeline with Stacks for Mysis species
#
# Steps:
#   1. Clone filtering
#   2. Adapter removal and trimming with cutadapt
#   3. De novo locus assembly with Stacks
#
# Assumptions:
#   - Raw files are named like: SAMPLE_R1_.fastq.gz and SAMPLE_R2_.fastq.gz
#   - clone_filter outputs paired files named like: SAMPLE.1.fq.gz / SAMPLE.2.fq.gz
#   - popmap file exists and is correctly formatted
#
# Adjust paths below before running.
###############################################################################

# ----------------------------- CONFIGURATION -------------------------------- #

RAW_DIR="$HOME/raw"
PROJECT_DIR="/scratch/project_2003448/full_dataset"

CLONE_DIR="$PROJECT_DIR/raw_clone_filter"
CUTADAPT_DIR="$PROJECT_DIR/cutadapt"
LOG_DIR="$PROJECT_DIR/LOGS"

STACKS_DIR="$PROJECT_DIR/stacks_denovo"
USTACKS_DIR="$STACKS_DIR/ustacks"
POPULATIONS_DIR="$STACKS_DIR/populations"

ADAPTERS_FILE="$CUTADAPT_DIR/adapters.fa"
POPMAP="$STACKS_DIR/popmap_mysis.txt"

THREADS=40
CUTADAPT_CORES=4

# Stacks parameters
USTACKS_M=3
USTACKS_MISMATCH=3
USTACKS_N=5
CSTACKS_N=5

# populations parameters
MIN_POPS=12
MIN_PROP=0.75

# ------------------------------ PREPARE DIRS -------------------------------- #

mkdir -p "$CLONE_DIR" \
         "$CUTADAPT_DIR" \
         "$LOG_DIR" \
         "$STACKS_DIR" \
         "$USTACKS_DIR" \
         "$POPULATIONS_DIR"

cd "$RAW_DIR"

###############################################################################
# 1. Clone filtering
###############################################################################

echo "Starting clone_filter..."

for r1 in *_R1_.fastq.gz; do
    sample="${r1%_R1_.fastq.gz}"
    r2="${sample}_R2_.fastq.gz"

    if [[ ! -f "$r2" ]]; then
        echo "WARNING: Missing mate pair for $r1 -> $r2 not found, skipping."
        continue
    fi

    echo "Processing clone_filter for sample: $sample"

    clone_filter \
        -1 "$r1" \
        -2 "$r2" \
        -i gzfastq \
        -o "$CLONE_DIR" \
        > "$CLONE_DIR/${sample}.log"
done

###############################################################################
# 2. Adapter removal and trimming with cutadapt
###############################################################################
# Assumes clone_filter output files are in CLONE_DIR and named:
#   sample.1.fq.gz
#   sample.2.fq.gz
###############################################################################

echo "Starting cutadapt..."

cd "$CLONE_DIR"

for r1 in *.1.fq.gz; do
    [[ -e "$r1" ]] || continue

    sample="${r1%.1.fq.gz}"
    r2="${sample}.2.fq.gz"

    if [[ ! -f "$r2" ]]; then
        echo "WARNING: Missing R2 file for $r1 -> $r2 not found, skipping."
        continue
    fi

    out_r1="$CUTADAPT_DIR/${sample}_trimmed.1.fq.gz"
    out_r2="$CUTADAPT_DIR/${sample}_trimmed.2.fq.gz"
    log_file="$LOG_DIR/${sample}_cutadapt.log"

    echo "Processing cutadapt for sample: $sample"

    cutadapt \
        -a "file:$ADAPTERS_FILE" \
        -A "file:$ADAPTERS_FILE" \
        -l 75 \
        -m 75 \
        --cores="$CUTADAPT_CORES" \
        -o "$out_r1" \
        -p "$out_r2" \
        "$r1" \
        "$r2" \
        > "$log_file"
done

###############################################################################
# 3. Stacks de novo assembly
###############################################################################
# Input to ustacks: trimmed R1 files
###############################################################################

echo "Starting ustacks..."

cd "$CUTADAPT_DIR"

id=1
for r1 in *_trimmed.1.fq.gz; do
    [[ -e "$r1" ]] || continue

    sample="${r1%_trimmed.1.fq.gz}"

    echo "Running ustacks for sample: $sample (ID: $id)"

    ustacks \
        -f "$r1" \
        -i "$id" \
        -t gzfastq \
        --name "$sample" \
        -m "$USTACKS_M" \
        -M "$USTACKS_MISMATCH" \
        -N "$USTACKS_N" \
        -p "$THREADS" \
        --force-diff-len \
        -o "$USTACKS_DIR"

    ((id+=1))
done

echo "Running cstacks..."
cstacks \
    -n "$CSTACKS_N" \
    -p "$THREADS" \
    -M "$POPMAP" \
    -P "$USTACKS_DIR"

echo "Running sstacks..."
sstacks \
    -P "$USTACKS_DIR" \
    -p "$THREADS" \
    -M "$POPMAP"

echo "Running tsv2bam..."
tsv2bam \
    -P "$USTACKS_DIR" \
    -t "$THREADS" \
    -M "$POPMAP" \
    -R "$CUTADAPT_DIR"

echo "Running gstacks..."
gstacks \
    -P "$USTACKS_DIR" \
    -t 8 \
    -M "$POPMAP" \
    --rm-pcr-duplicates

echo "Running populations..."
populations \
    -P "$USTACKS_DIR" \
    -O "$POPULATIONS_DIR" \
    --popmap "$POPMAP" \
    -p "$MIN_POPS" \
    -r "$MIN_PROP" \
    --write-random-snp \
    --hwe \
    --fstats \
    -k \
    --structure \
    --genepop \
    --vcf \
    --plink \
    --phylip \
    --treemix \
    --fasta-loci \
    --fasta-samples-raw \
    --phylip-var-all \
    --phylip-var

echo "Pipeline finished successfully."
