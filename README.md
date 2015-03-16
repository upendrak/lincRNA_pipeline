# lincRNA_pipeline
This pipeline generates long non-coding RNA from TopHat-Cufflinks-Cuffcompare output

#### Usage Summary

The following can be accessed by running `sh lincRNA_pipeline.sh -h`:

`USAGE : sh lincRNA_pipeline.sh` 

              `-c    /path/to/cuffcompare_output file`
              `-g    /path/to/reference genome file`
              `-r    /path/to/reference CDS file`
              `-b    /path/to/RNA file`

#### Example

`sh lincRNA_pipeline.sh -c cuffcompare_out_annot_no_annot.combined.gtf -g Brapa_sequence_v1.2_genome.fa -r Brassica_rapa_v1.2_cds.fa -b TE_RNA_transcripts.fa`

