
count_reads() {
    for fq in "$@"; do
        if [[ "$fq" == *.gz ]]; then
            echo $(( $(zcat "$fq" | wc -l) / 4 ))
        else
            echo $(( $(wc -l < "$fq") / 4 ))
        fi
    done
}

pair_fq_count_check () {
    fastq1="$1"
    fastq2="$2"
    output_file="$3"

    reads1=$(count_reads "$fastq1")
    reads2=$(count_reads "$fastq2")
    echo $fastq1 $reads1
    echo $fastq2 $reads2

    if [[ "$reads1" -eq "$reads2" ]]; then
        echo "$reads1" > "$output_file"
    fi
}
# Usage:
#   fastq_len_qc <fastq|fastq.gz> <out_prefix> [binsize] [break_min] [break_max]
# Examples:
#   fastq_len_qc reads.fastq.gz out/readlen
#   fastq_len_qc reads.fastq.gz out/readlen 25
#   fastq_len_qc reads.fastq.gz out/readlen 25 10 500   # clamp <10 to 10, >500 to 500 for histogram

fastq_len_qc() {
  local fq="$1"
  local prefix="$2"
  local binsize="${3:-50}"
  local break_min="${4:-}"
  local break_max="${5:-}"

  if [[ -z "$fq" || -z "$prefix" ]]; then
    echo "Usage: fastq_len_qc <fastq|fastq.gz> <out_prefix> [binsize] [break_min] [break_max]" >&2
    return 2
  fi
  if [[ ! -f "$fq" ]]; then
    echo "ERROR: file not found: $fq" >&2
    return 2
  fi
  command -v seqkit >/dev/null 2>&1 || { echo "ERROR: seqkit not found in PATH" >&2; return 127; }
  command -v Rscript >/dev/null 2>&1 || { echo "ERROR: Rscript not found in PATH" >&2; return 127; }

  # validate break args (either both empty OR both provided)
  if [[ -n "$break_min" || -n "$break_max" ]]; then
    if [[ -z "$break_min" || -z "$break_max" ]]; then
      echo "ERROR: break_min and break_max must be provided together, or omit both." >&2
      return 2
    fi
  fi

  local stats_out="${prefix}.seqkit_stats_a.txt"
  local len_out="${prefix}.fx2tab_l.tsv"   # keep ALL columns from fx2tab -l
  local q_out="${prefix}.read_length_quantiles.tsv"
  local hist_png="${prefix}.read_length_hist_bin${binsize}.png"

  # If breaks provided, reflect in output filename (so replot doesn't overwrite)
  if [[ -n "$break_min" && -n "$break_max" ]]; then
    hist_png="${prefix}.read_length_hist_bin${binsize}.clamp${break_min}-${break_max}.png"
  fi

  echo "==> [1/3] seqkit stats -a (write + print): $stats_out" >&2
  if [[ -s "$stats_out" ]]; then
    echo "    Found existing $stats_out (non-empty). Skipping regeneration." >&2
  else
    # print AND write
    seqkit stats -a "$fq" | tee "$stats_out"
  fi

  echo "==> [2/3] seqkit fx2tab -l (keep full output): $len_out" >&2
  if [[ -s "$len_out" ]]; then
    echo "    Found existing $len_out (non-empty). Skipping regeneration." >&2
  else
    seqkit fx2tab -l "$fq" > "$len_out"
  fi

  echo "==> [3/3] Rscript: histogram + quantiles (binsize=$binsize${break_min:+, break_min=$break_min}${break_max:+, break_max=$break_max})" >&2
  Rscript --vanilla - "$len_out" "$prefix" "$binsize" "${break_min:-NA}" "${break_max:-NA}" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
fx2tab_file <- args[1]
prefix      <- args[2]
binsize     <- as.numeric(args[3])
break_min_s <- args[4]
break_max_s <- args[5]

if (is.na(binsize) || binsize <= 0) stop("binsize must be a positive number")

use_clamp <- !(toupper(break_min_s) == "NA" || toupper(break_max_s) == "NA")
break_min <- suppressWarnings(as.numeric(break_min_s))
break_max <- suppressWarnings(as.numeric(break_max_s))

if (use_clamp) {
  if (is.na(break_min) || is.na(break_max)) stop("break_min/break_max must be numeric when provided")
  if (break_min >= break_max) stop("break_min must be < break_max")
}

# Fast extract last column (lengths)
cmd <- sprintf("awk -F'\\t' '{print $NF}' %s", shQuote(fx2tab_file))
lens <- suppressWarnings(as.integer(system(cmd, intern=TRUE)))
lens <- lens[!is.na(lens)]
if (length(lens) == 0) stop("No lengths parsed from last column of: ", fx2tab_file)

# Quantiles (based on RAW lengths, not clamped)
probs <- c(0, 0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 1.0)
qs <- as.numeric(quantile(lens, probs=probs, names=FALSE, type=7))

q_out <- paste0(prefix, ".read_length_quantiles.tsv")
q_df <- data.frame(prob=probs, quantile=qs)
write.table(q_df, file=q_out, sep="\t", quote=FALSE, row.names=FALSE)

cat("\n# Read length quantiles (raw lengths)\n")
print(q_df, row.names=FALSE)

# Histogram (optionally clamp extremes into edge buckets)
lens_hist <- lens
if (use_clamp) {
  lens_hist[lens_hist < break_min] <- break_min
  lens_hist[lens_hist > break_max] <- break_max

  # breaks must include endpoints; avoid dropping values
  # Use: [break_min, break_min+binsize, ..., break_max]
  br <- seq(break_min, break_max, by=binsize)
  if (tail(br, 1) != break_max) br <- c(br, break_max)
  breaks <- unique(br)

  png_out <- paste0(prefix, ".read_length_hist_bin", binsize, ".clamp", break_min, "-", break_max, ".png")
  main_t <- paste0("Read length distribution (bin=", binsize, ", clamp <", break_min, " and >", break_max, ")")
} else {
  maxx <- max(lens_hist)
  breaks <- seq(0, maxx + binsize, by=binsize)

  png_out <- paste0(prefix, ".read_length_hist_bin", binsize, ".png")
  main_t <- paste0("Read length distribution (bin=", binsize, ")")
}

png(png_out, width=1200, height=800)
hist(lens_hist, breaks=breaks,
     right=FALSE, include.lowest=TRUE,
     main=main_t,
     xlab="Read length", ylab="Count")
dev.off()

cat("\n# Wrote:\n")
cat("  ", q_out, "\n", sep="")
cat("  ", png_out, "\n", sep="")
RSCRIPT

  echo "==> Done. Outputs:" >&2
  echo "  $stats_out" >&2
  echo "  $len_out" >&2
  echo "  $q_out" >&2
  echo "  $hist_png" >&2
}
