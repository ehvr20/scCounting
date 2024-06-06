#!/bin/bash

set -euo pipefail

# Configuration
CONFIG=../configuration.json
OUTDIR=../results
TEMPDIR=../temp

# Parse command-line options
# -i whether intronic reads should be counted
intron_flag=""
while getopts "i" opt; do
  case ${opt} in
    i )
      intron_flag='--soloFeatures GeneFull_Ex50pAS'
      ;;
    \? )
      echo "Usage: $0 [-i] manifest"
      exit 1
      ;;
  esac
done
shift $((OPTIND - 1))

# Check for manifest argument
if [ $# -eq 0 ]; then
  echo "No manifest given." >&2
  exit 1
else
  manifest=$1
fi

# Extract sample name
sample=$(basename "$manifest" | cut -d '_' -f 1)

# Make output directory
resultdir="$OUTDIR/outs/$sample"
mkdir -p "$resultdir"

# Read configuration values
star_image=$(jq -r ".images.star" "$CONFIG")
threads=$(jq -r ".resources.cores" "$CONFIG")
genome_dir=$(jq -r ".reference_genome" "$CONFIG")
gtf_file=$(jq -r ".gtf_file" "$CONFIG")

# Run STAR alignment
singularity run "$star_image" STAR \
  --runMode alignReads \
  --runThreadN "$threads" \
  --quantMode GeneCounts \
  --soloType SmartSeq --soloStrand Unstranded --soloUMIdedup NoDedup \
  $intron_flag \
  --outFilterScoreMin 30 \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix "$resultdir/" \
  --readFilesCommand zcat \
  --genomeDir "$genome_dir" \
  --sjdbGTFfile "$gtf_file" \
  --readFilesManifest "$(realpath "$manifest")" || exit 1

# Create symbolic links for outputs
ln -sr "$resultdir/Solo.out/*Gene*/raw" "$OUTDIR/${sample}_raw"
ln -sr "$resultdir/Solo.out/*Gene*/filtered" "$OUTDIR/${sample}_filtered"

