**NAME**

lincRNA_pipeline

**DESCRIPTION**

This document describes a pipeline that integrates several Python script and linux commands to generate long non-coding RNA (lncRNA) from TopHat-Cufflinks-Cuffcompare output

**USUAGE NOTES**

The pipeline can be run by using the script `lincRNA_pipeline.sh`

./lincRNA_pipeline.sh 

              	-c    /path/to/cuffcompare_output file
                -g    /path/to/reference genome file
                -r    /path/to/reference CDS file
                -b    /path/to/RNA file

**EXAMPLE**

./lincRNA_pipeline.sh -c cuffcompare_out_annot_no_annot.combined.gtf -g Brapa_sequence_v1.2_genome.fa -r Brassica_rapa_v1.2_cds.fa -b TE_RNA_transcripts.fa

