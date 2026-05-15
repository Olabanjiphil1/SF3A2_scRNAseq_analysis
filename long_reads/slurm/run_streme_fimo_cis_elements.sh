#!/bin/bash
#SBATCH --job-name=cis-element-finder
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=cis-element-finder.out
#SBATCH --error=cis-element-finder.err

set -euo pipefail

# Edit these paths or export them before running.
GENOME="${GENOME:-data/long_reads/genome/PlasmoDB-68_Pfalciparum3D7_Genome.fasta}"
OUT="${OUT:-results/long_reads/motif}"

# Add MEME Suite to PATH if needed:
# export PATH=/path/to/meme/bin:/path/to/meme/libexec/meme-5.5.9:$PATH

for tool in bedtools streme fimo; do
  command -v "$tool" >/dev/null 2>&1 || {
    echo "ERROR: $tool not found in PATH" >&2
    exit 1
  }
done

mkdir -p "$OUT"

echo "Genome: $GENOME"
echo "Output: $OUT"

for window in donor acceptor_bpaware branchpoint_focus; do
  bedtools getfasta -fi "$GENOME" -bed "$OUT/SF3A2_positive_${window}.bed" -s -name -fo "$OUT/SF3A2_positive_${window}.fa"
  bedtools getfasta -fi "$GENOME" -bed "$OUT/NF54_background_${window}.bed" -s -name -fo "$OUT/NF54_background_${window}.fa"

  streme \
    --p "$OUT/SF3A2_positive_${window}.fa" \
    --n "$OUT/NF54_background_${window}.fa" \
    --dna \
    --oc "$OUT/STREME_${window}"

  fimo \
    --oc "$OUT/FIMO_${window}_positive" \
    "$OUT/STREME_${window}/streme.xml" \
    "$OUT/SF3A2_positive_${window}.fa"

  fimo \
    --oc "$OUT/FIMO_${window}_background" \
    "$OUT/STREME_${window}/streme.xml" \
    "$OUT/NF54_background_${window}.fa"
done
