#!/bin/bash

> ../temp/launch_commands.csv

dataset_manifest_file=$(jq -r ".manifest" ../configuration.json)
export dataset_manifest_file

all_files=$(find ../data/ -iname "*.fastq.gz" | grep '_S[0-9]*_L[0-9]*_[R,I][0-9]_[0-9]*')
export all_files

samples=$(echo "$all_files" | xargs -n 1 -I {} basename {} | awk -F'_S[0-9]_' '{print $1}' | cut -d'_' -f 1 | sort -u)

function get_count_command () {

	sample="$1"

	if [ -L "../results/${sample}_raw" ] && [ -e "../results/${sample}_raw" ]; then
		# echo "Skipping $sample..."
		return
	fi

	sample_files=$(echo "$all_files" | grep "$sample")

	refid=$(echo "$sample" | grep -o "v[0-9]-[0-9]*")

	suspension=$(jq -r ".\"$refid\".Suspension" "$dataset_manifest_file")
	intron_flag=""
	if $(echo $suspension | grep -qi 'sn'); then
		intron_flag='-i'
	fi

	technology=$(jq -r ".\"$refid\".Technology" "$dataset_manifest_file")

	# Construct a manifest and return the line to be run to run the counts
	case "${technology,,}" in
		*drop*)
			temp_manifest=../temp/${sample}_manifest.tsv

			read1=$(echo "$sample_files" | grep "_R1_")
			read2=$(echo "$sample_files" | grep "_R2_")
			cell_ids=$(echo "$read2" | awk -F'_S[0-9]_' '{print $1}')

			cbumi_ratio=$(jq -r ".\"$refid\".\"BC+UMI\"" "$dataset_manifest_file")
			barcode_flags=""
			if [ ! $cbumi_ratio = "null" ]; then
				cb_len=$(echo "$cbumi_ratio" | cut -d '+' -f 1)
				umi_len=$(echo "$cbumi_ratio" | cut -d '+' -f 2)
				barcode_flags="-c $cb_len -u $umi_len"
			fi

			paste -d '\t' <(echo "${read2}") <(echo "${read1}") <(echo "${cell_ids}") > $temp_manifest

			echo "utils/count_dropseq_reads.sh $intron_flag $barcode_flags $temp_manifest" >> ../temp/launch_commands.csv
			;;

		*smart*)
			temp_manifest=../temp/${sample}_manifest.tsv

			read1=$(echo "$sample_files" | grep "_R1_")
			read2=$(echo "$sample_files" | grep "_R2_")
			cell_ids=$(echo "$read2" | awk -F'_S[0-9]_' '{print $1}')

			# If no read 1 given, replace with a dash
			if [ -z "$read1"]; then
				read1=$(yes "-" | head -n $(echo "$read2" | wc -l))
			fi

			# Output all files into tsv manifest of 3 columns: 	read2	read1	sample_id
			paste -d '\t' <(echo "${read2}") <(echo "${read1}") <(echo "${cell_ids}") > $temp_manifest

			echo "utils/count_smartseq_reads.sh $intron_flag $temp_manifest" >> ../temp/launch_commands.csv
			;;

		*10x*)
			temp_manifest=../temp/${sample}_manifest.csv

			tenx_chemistry=$(jq -r ".\"$refid\".\"10x_chemistry\"" "$dataset_manifest_file")
			chemistry_arg=""
                        if [ ! ${tenx_chemistry} = "null" ]; then
                                chemistry_arg="-c $tenx_chemistry"
                        fi

			echo "$sample_files" | xargs -n 1 -I {} basename {} | awk -F'_S[0-9]_' '{print $1}' | sort -u | paste -sd "," - > $temp_manifest

			echo "utils/count_tenx_reads.sh $intron_flag $chemistry_arg $temp_manifest" >> ../temp/launch_commands.csv
			;;
	esac

}
export -f get_count_command

echo "$samples" | xargs -n 1 -P 4 -I {} bash -c "get_count_command {}"

# while read command; do
# 	echo $command
# done < ../temp/launch_commands.csv
# exit 0

SLURM_JOB_NAME="count_reads"
SLURM_ACCOUNT=$(jq -r ".slurm.cpu.account" ../configuration.json)
SLURM_PARTITION=$(jq -r ".slurm.cpu.partition" ../configuration.json)
SLURM_LOGFILE="../logs/${SLURM_JOB_NAME}-%A_%a.log"
SLURM_CORES=$(jq -r ".resources.cores" ../configuration.json)

function job_script () {
	command=$(head -n $SLURM_ARRAY_TASK_ID "$1" | tail -n 1)
	bash $command
}
export -f job_script

array_len=$(wc -l < ../temp/launch_commands.csv)

if [ $array_len -eq 0 ]; then
	echo "All samples already processed"
	exit 0
else
	sbatch \
	-J $SLURM_JOB_NAME \
	-A $SLURM_ACCOUNT \
	-p $SLURM_PARTITION \
	-o $SLURM_LOGFILE \
	-t 24:00:00 \
	-c $SLURM_CORES \
	-a 1-$array_len \
	--wrap="job_script '../temp/launch_commands.csv'"
fi

