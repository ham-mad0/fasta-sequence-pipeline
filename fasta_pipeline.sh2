#!/usr/bin/env bash
# =====================================================================
#  FASTA FAST-TRACK PIPELINE  (A to Z)
#  Download -> Stats -> GC content -> Window GC -> N50 -> Reverse
#  complement -> (optional) BWA alignment -> BAM -> VCF
#
#  Usage:
#    bash fasta_pipeline.sh NC_000913.3                      (NCBI se download)
#    bash fasta_pipeline.sh my_genome.fasta                  (apni FASTA file)
#    bash fasta_pipeline.sh NC_002695.2 NC_000913.3.fasta    (reference ke saath)
#
#  Window size badalna ho:  WINDOW=5000 bash fasta_pipeline.sh NC_000913.3
# =====================================================================
set -euo pipefail

INPUT="${1:-}"
REF="${2:-}"
WINDOW="${WINDOW:-10000}"

if [ -z "$INPUT" ]; then
    echo "Ghalti! Sahi tarika:"
    echo "  bash fasta_pipeline.sh <ACCESSION ya fasta_file> [reference.fasta]"
    exit 1
fi

# ---------------------------------------------------------------------
# 1. FASTA hasil karna (file ho to wohi, warna NCBI se download)
# ---------------------------------------------------------------------
echo "=== [1/7] FASTA hasil ho rahi hai ==="
if [ -f "$INPUT" ]; then
    FASTA="$INPUT"
    NAME="$(basename "$INPUT")"
    NAME="${NAME%.*}"
    echo "Local file istemal ho rahi hai: $FASTA"
else
    NAME="$INPUT"
    FASTA="${NAME}.fasta"
    echo "NCBI se download: $NAME"
    curl -fSL --retry 3 \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${NAME}&rettype=fasta&retmode=text" \
      -o "$FASTA" || { echo "Download fail. Internet ya accession check karein."; exit 1; }
fi

# Check: waqai FASTA hai?
if ! head -c 1 "$FASTA" | grep -q '>'; then
    echo "Ghalti: yeh FASTA file nahi lagti. Pehli line dekhein: $(head -n 1 "$FASTA")"
    exit 1
fi

OUT="results_${NAME}"
mkdir -p "$OUT"

# ---------------------------------------------------------------------
# 2. Har sequence ke statistics (length, A/T/G/C/N, GC%)
# ---------------------------------------------------------------------
echo "=== [2/7] Per-sequence statistics ==="
printf "sequence_id\tlength\tA\tT\tG\tC\tN\tGC_percent\n" > "$OUT/per_sequence_stats.tsv"
awk '
function flush(   tot, gcp) {
    if (id == "") return
    tot = a + t + g + c
    gcp = (tot > 0) ? (g + c) * 100 / tot : 0
    printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n", id, len, a, t, g, c, n, gcp
}
/^>/ { flush(); id = substr($1, 2); len = a = t = g = c = n = 0; next }
{
    s = toupper($0); gsub(/[^A-Z]/, "", s); len += length(s)
    x = s; a += gsub(/A/, "", x)
    x = s; t += gsub(/T/, "", x)
    x = s; g += gsub(/G/, "", x)
    x = s; c += gsub(/C/, "", x)
    x = s; n += gsub(/N/, "", x)
}
END { flush() }
' "$FASTA" >> "$OUT/per_sequence_stats.tsv"

# ---------------------------------------------------------------------
# 3. Kul summary (sequences, length, GC%, N50)
# ---------------------------------------------------------------------
echo "=== [3/7] Overall summary (GC content, N50) ==="
tail -n +2 "$OUT/per_sequence_stats.tsv" | cut -f2 | sort -rn > "$OUT/.lengths.tmp"

tail -n +2 "$OUT/per_sequence_stats.tsv" | awk -F'\t' '
{ n++; len += $2; A += $3; T += $4; G += $5; C += $6; N += $7
  if (min == "" || $2 < min) min = $2
  if ($2 > max) max = $2 }
END {
    tot = A + T + G + C
    gc = (tot > 0) ? (G + C) * 100 / tot : 0
    at = (tot > 0) ? (A + T) * 100 / tot : 0
    printf "total_sequences\t%d\ntotal_length_bp\t%d\nlongest_bp\t%d\nshortest_bp\t%d\n", n, len, max, min
    printf "A\t%d\nT\t%d\nG\t%d\nC\t%d\nN\t%d\n", A, T, G, C, N
    printf "GC_percent\t%.2f\nAT_percent\t%.2f\n", gc, at
}' > "$OUT/summary.tsv"

# N50
awk -v total="$(awk '{s+=$1} END{print s+0}' "$OUT/.lengths.tmp")" '
{ cum += $1; if (cum >= total / 2) { print "N50_bp\t" $1; exit } }
' "$OUT/.lengths.tmp" >> "$OUT/summary.tsv"
rm -f "$OUT/.lengths.tmp"

cat "$OUT/summary.tsv"

# ---------------------------------------------------------------------
# 4. Sliding-window GC content (plot ya GC skew dekhne ke liye)
# ---------------------------------------------------------------------
echo "=== [4/7] Window GC content (window = ${WINDOW} bp) ==="
awk -v W="$WINDOW" '
function flush_rec(   i, w, x, gcn, atn, tot) {
    if (id == "") return
    for (i = 1; i <= length(seq); i += W) {
        w = toupper(substr(seq, i, W))
        x = w; gcn = gsub(/[GC]/, "", x)
        x = w; atn = gsub(/[AT]/, "", x)
        tot = gcn + atn
        printf "%s\t%d\t%d\t%.2f\n", id, i, i + length(w) - 1, (tot > 0) ? gcn * 100 / tot : 0
    }
}
BEGIN { print "sequence_id\tstart\tend\tGC_percent" }
/^>/ { flush_rec(); id = substr($1, 2); seq = ""; next }
{ seq = seq $0 }
END { flush_rec() }
' "$FASTA" > "$OUT/gc_windows.tsv"
echo "Saved: $OUT/gc_windows.tsv"

# ---------------------------------------------------------------------
# 5. Reverse complement
# ---------------------------------------------------------------------
echo "=== [5/7] Reverse complement ==="
awk '
/^>/ { if (h != "") print h "\t" s; h = $0; s = ""; next }
{ s = s $0 }
END { if (h != "") print h "\t" s }
' "$FASTA" | while IFS=$'\t' read -r header seq; do
    echo "${header} [reverse_complement]"
    echo "$seq" | rev | tr 'ACGTacgt' 'TGCAtgca' | fold -w 70
done > "$OUT/${NAME}_revcomp.fasta"
echo "Saved: $OUT/${NAME}_revcomp.fasta"

# ---------------------------------------------------------------------
# 6 & 7. Optional: reference ke saath alignment (BWA -> BAM -> VCF)
# ---------------------------------------------------------------------
if [ -n "$REF" ]; then
    if [ ! -f "$REF" ]; then
        echo "Ghalti: reference file nahi mili: $REF"
        exit 1
    fi
    if command -v bwa >/dev/null 2>&1 && command -v samtools >/dev/null 2>&1; then
        echo "=== [6/7] BWA alignment -> sorted BAM ==="
        [ -f "${REF}.bwt" ] || bwa index "$REF"
        bwa mem "$REF" "$FASTA" | samtools sort -o "$OUT/${NAME}_sorted.bam" -
        samtools index "$OUT/${NAME}_sorted.bam"
        samtools flagstat "$OUT/${NAME}_sorted.bam" > "$OUT/${NAME}_flagstat.txt"
        cat "$OUT/${NAME}_flagstat.txt"

        if command -v bcftools >/dev/null 2>&1; then
            echo "=== [7/7] Variant calling -> VCF ==="
            [ -f "${REF}.fai" ] || samtools faidx "$REF"
            bcftools mpileup -f "$REF" "$OUT/${NAME}_sorted.bam" \
              | bcftools call -mv -Ov -o "$OUT/${NAME}.vcf"
            echo "Variants mile: $(grep -vc '^#' "$OUT/${NAME}.vcf")"
        else
            echo "bcftools install nahi hai. Install: mamba install -n wgs -c bioconda bcftools"
        fi
    else
        echo "bwa/samtools nahi mile. Pehle chalayein: mamba activate wgs"
    fi
else
    echo "=== [6-7/7] Reference nahi di, alignment skip ==="
fi

echo ""
echo "=== KHATAM! Sab nateeje yahan hain: $OUT/ ==="
ls -1 "$OUT"
