#! /bin/bash

nextflow run nf-core/ampliseq -r 2.8.0 --input Samplesheet.tsv --metadata PAN_Metadata.tsv --outdir PAN_ampliseq -profile singularity --picrust --metadata_category Sex,Injury,Treatment,Timepoint --FW_primer AGAGTTTGATYMTGGCTCAG --RV_primer ATTACCGCGGCKGCTGG --trunclenf 275 --trunclenr 265 --skip_cutadapt --max_cpus 12 -resume