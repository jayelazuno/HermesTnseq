

#!/bin/bash
## Title: Function to concatenate multiplexed FASTQ files for a single sample
## Author: Joshua Ayelazuno
## Date: 08-07-2025
#
# Go to the FASTQ file directory
cd /Users/jayelazuno/workspace/HermesTnseq/data || { echo "Directory not found!"; exit 1; }

concat_sample_reads() {
  sample_prefix=$1              # e.g., S1
  output_file="${sample_prefix}_R1.fastq.gz"
  input_files=$(ls ${sample_prefix}[a-e]_R1.fastq.gz 2>/dev/null)

  if [ -z "$input_files" ]; then
    echo " No files found for sample ${sample_prefix}"
    return 1
  fi

  echo " Merging files for ${sample_prefix}..."
  echo "$input_files"
  cat $input_files > $output_file
  echo " Output written to: $output_file"
  echo ""
}

# Loop over all sample groups: S1 through S8
for i in {1..8}; do
  concat_sample_reads "S${i}"
done

echo " All samples merged."
