#!/bin/bash
	
	for file in $(<filelist_PE_QC.txt)
do
	java -jar $TRIM/trimmomatic  PE -threads 4 -trimlog ${file}_QC_PE_log.fasta ${file}_L001_R1_001.fastq ${file}_L001_R2_001.fastq ${file}_L001_R1_001_pTrim.fastq ${file}_L001_R1_001_uTrim.fastq ${file}_L001_R2_001_pTrim.fastq ${file}_L001_R2_001_uTrim.fastq ILLUMINACLIP:AdapterRef16S.fa:2:30:12 LEADING:20 TRAILING:20 MINLEN:125

done


