#!/bin/bash

set -euo pipefail

# Configuration
CONFIG=../configuration.json
OUTDIR=../results
TEMPDIR=../temp

# Parse command-line options
# -c length of cell barcode
# -u length of umi
# -i whether intronic reads should be counted
intron_flag=""
cb_len="12"
umi_len="8"
while getopts "ic:u:" opt; do
  case ${opt} in
    i )
      intron_flag='--soloFeatures GeneFull_Ex50pAS'
      ;;
    c )
      cb_len=$OPTARG
      ;;
    u )
      umi_len=$OPTARG
      ;;
    \? )
      echo "Usage: $0 [-i] [-c 12] [-u 8] manifest"
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

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

# Define cb_start and umi_start
cb_start=1
umi_start=$((cb_start + cb_len))

# Run STAR alignment
singularity run $star_image STAR \
    --runMode=alignReads \
    --soloType CB_UMI_Simple --soloStrand Forward \
    --soloCBstart "$cb_start" --soloCBlen "$cb_len" --soloUMIstart "$umi_start" --soloUMIlen "$umi_len" --soloCBwhitelist None \
    --soloBarcodeReadLength 0 \
    --soloCellFilter EmptyDrops_CR \
    --soloUMIdedup 1MM_CR \
    $intron_flag \
    --quantMode GeneCounts \
    --runThreadN $threads \
    --limitOutSJcollapsed 20000000 \
    --outSAMtype BAM Unsorted \
    --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
    --outFileNamePrefix "$resultdir/" \
    --readFilesCommand zcat \
    --genomeDir "$genome_dir" \
    --sjdbGTFfile "$gtf_file" \
    --readFilesManifest $(realpath "$manifest") || exit 1

# Create symbolic links for outputs
ln -sr "$resultdir/Solo.out/*Gene*/raw" "$OUTDIR/${sample}_raw"
ln -sr "$resultdir/Solo.out/*Gene*/filtered" "$OUTDIR/${sample}_filtered"

