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
intron_flag="--include-introns=false"
chemistry="--chemistry=auto"
while getopts "ic:" opt; do
  case ${opt} in
    i )
      intron_flag='--include-introns=true'
      ;;
    c )
      chemistry="--chemistry=${OPTARG}"
      ;;
    \? )
      echo "Usage: $0 [-i] [-c] manifest"
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
output_directory="$OUTDIR/outs"
mkdir -p "$output_directory"

resultdir="$OUTDIR/outs/$sample"

# Read configuration values
cellranger_image=$(jq -r ".images.cellranger" "$CONFIG")
memory=$(jq -r ".resources.memory" "$CONFIG")
transcriptome=$(jq -r ".cr_transcriptome" "$CONFIG")

workdir=$(pwd)

sample_ids=$(cat $manifest)
fastqs_dir=$(realpath ../data)

cd "$output_directory"
# Run STAR alignment
singularity run "${cellranger_image}" \
 cellranger count --transcriptome $transcriptome --fastqs $fastqs_dir \
  --sample "$sample_ids" --id $sample \
  ${intron_flag} \
  ${multiomic_flag} \
  --localmem $memory || exit 1

cd "$workdir"
# Create symbolic links for outputs
ln -sr "$resultdir/outs/raw_feature_bc_matrix" "$OUTDIR/${sample}_raw"
ln -sr "$resultdir/outs/filtered_feature_bc_matrix" "$OUTDIR/${sample}_filtered"

rm -r "$resultdir/SC_RNA_COUNTER_CS"
