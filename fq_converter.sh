#!/bin/bash

# This will change all the *.fq files to .fastq
for file in *.fq; do
    # Check if the corresponding .fastq file already exists
    if [ ! -f "${file%.fq}.fastq" ]; then
        mv "$file" "${file%.fq}.fastq"
    else
        echo "Skipped: ${file%.fq}.fastq already exists."
    fi
done

# Now gzip the .fastq files, if already done then will skip this time consuming step
for file in *.fastq *.fastq.gz; do
  if [[ $file == *.fastq ]]; then
    # Compress .fastq files
    pigz -p10 "$file"
  elif [[ $file == *.fastq.gz ]]; then
    # Skip .fastq.gz files
    echo "Skipping $file as it is already compressed."
  fi
done

# Now we want to convert the *.1.fastq.gz to *_R1_001.fastq.gz
for file in *.fastq.gz; do
  if [[ $file == *.1.fastq.gz ]]; then
    base=${file%.1.fastq.gz}
    newname="${base}_R1_001.fastq.gz"
  elif [[ $file == *.2.fastq.gz ]]; then
    base=${file%.2.fastq.gz}
    newname="${base}_R2_001.fastq.gz"
  else
    continue
  fi
  mv "$file" "$newname"
done


