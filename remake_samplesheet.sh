#!/bin/bash

# Create or overwrite the Samplesheet.tsv file with headers
echo -e "sampleID\tforwardReads\treverseReads" > Samplesheet.tsv

# Loop through all forward read files
for forward in *_R1_*.fastq.gz; do
    # Extract the sample ID
    sampleID=$(echo $forward | sed 's/_R1_.*//')

    # Construct the reverse read filename
    reverse="${sampleID}_R2_001.fastq.gz"

    # Check if the reverse read file exists
    if [ -f "$reverse" ]; then
        # Add the entry to the Samplesheet
        echo -e "${sampleID}\t$(pwd)/$forward\t$(pwd)/$reverse" >> Samplesheet.tsv
    else
        echo "Reverse read file for $sampleID not found."
    fi
done

