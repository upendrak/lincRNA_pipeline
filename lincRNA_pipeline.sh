#!/bin/bash

# Script to process cuffcompare output file to generate lincRNA
# Usage: 
# sh lincRNA_pipeline.sh -c cuffcompare_out_annot_no_annot.combined.gtf -g Brapa_sequence_v1.2.fa -r Brassica_rapa_v1.2.cds -b TE_RNA_transcripts.fa

while getopts ":b:c:g:hr:" opt; do
  case $opt in
    b)
      blastfile=$OPTARG
      ;;
    c)
      comparefile=$OPTARG
      ;;
    h)	echo 	"USAGE : sh lincRNA_pipeline.sh 
		          -c 	</path/to/cuffcompare_output file>
	 		  -g 	</path/to/reference genome file>
                          -r    </path/to/reference CDS file>
			  -b 	</path/to/RNA file>"
      exit 1
      ;;
    g)
      referencegenome=$OPTARG
      ;;
    r)
     referenceCDS=$OPTARG
      ;;	
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

wget -O- https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz | tar xzvf -

makeblastdb -in $blastfile -dbtype nucl -out $blastfile.blast.out

makeblastdb -in $referenceCDS -dbtype nucl -out $referenceCDS.blast.out

grep '"u"' $comparefile | \
	gffread -w transcripts_u.fa -g $referencegenome - && \
	python2.7 get_gene_length_filter.py transcripts_u.fa transcripts_u_filter.fa && \
	TransDecoder-2.0.1/TransDecoder.LongOrfs -t transcripts_u_filter.fa

sed 's/ .*//' transcripts_u_filter.fa | sed -ne 's/>//p' \
	> transcripts_u_filter.fa.genes

cd transcripts_u_filter.fa.transdecoder_dir

sed 's/|.*//' longest_orfs.cds | sed -ne 's/>//p' | uniq \
	> longest_orfs.cds.genes

grep -v -f longest_orfs.cds.genes ../transcripts_u_filter.fa.genes \
	> longest_orfs.cds.genes.not.genes 


sed 's/^/>/' longest_orfs.cds.genes.not.genes \
	> temp && mv temp longest_orfs.cds.genes.not.genes

python ../extract_sequences.py longest_orfs.cds.genes.not.genes ../transcripts_u_filter.fa longest_orfs.cds.genes.not.genes.fa

sed 's/ /./' longest_orfs.cds.genes.not.genes.fa \
  > temp && mv temp longest_orfs.cds.genes.not.genes.fa

blastn -query longest_orfs.cds.genes.not.genes.fa -db ../$blastfile.blast.out -out longest_orfs.cds.genes.not.genes.fa.blast.out -outfmt 6

python ../filter_sequences.py longest_orfs.cds.genes.not.genes.fa.blast.out longest_orfs.cds.genes.not.genes.fa.blast.out.filtered

grep ">" longest_orfs.cds.genes.not.genes.fa | sed 's/>//' \
  > longest_orfs.cds.genes.not.genes_only

python ../fasta_remove.py longest_orfs.cds.genes.not.genes.fa.blast.out.filtered longest_orfs.cds.genes.not.genes_only lincRNA.genes

sed 's/^/>/' lincRNA.genes > temp && mv temp lincRNA.genes

python ../extract_sequences-1.py lincRNA.genes longest_orfs.cds.genes.not.genes.fa lincRNA.genes.fa

cap3 lincRNA.genes.fa

cat lincRNA.genes.fa.cap.singlets lincRNA.genes.fa.cap.contigs \
  > lincRNA_genes_non_redundant.fa

grep ">" lincRNA_genes_non_redundant.fa \
  > lincRNA_gene_non_redundant.genes_only


blastn -query lincRNA_genes_non_redundant.fa -db ../$referenceCDS.blast.out -out lincRNA_non_redundant.fa_cds_blast.out -evalue 1e-30 -outfmt 6

python ../linc_RNA_filter.py lincRNA_non_redundant.fa_cds_blast.out lincRNA_non_redundant.fa_cds_blast.out.filtered

grep -v -f lincRNA_non_redundant.fa_cds_blast.out.filtered lincRNA_gene_non_redundant.genes_only \
  > lincRNA_non_redundant_filtered.genes

python ../extract_sequences-1.py lincRNA_non_redundant_filtered.genes lincRNA_genes_non_redundant.fa lincRNA_non_redundant_filtered.genes.fa

grep "^>" lincRNA_non_redundant_filtered.genes.fa \
  > lincRNA_non_redundant_filtered.genes.only

makeblastdb -in lincRNA_non_redundant_filtered.genes.fa -dbtype nucl -out lincRNA_non_redundant_filtered_blast.out

blastn -query lincRNA_non_redundant_filtered.genes.fa -db lincRNA_non_redundant_filtered_blast.out -out lincRNA_non_redundant_filtered_blast.out.blast.out -evalue 1e-30 -outfmt 6

python ../linc_RNA_filter-1.py lincRNA_non_redundant_filtered_blast.out.blast.out lincRNA_non_redundant_filtered_blast.out.blast.out.uniq

grep -v -f lincRNA_non_redundant_filtered_blast.out.blast.out.uniq lincRNA_non_redundant_filtered.genes.only \
  > lincRNA_non_redundant_filtered.genes.only_filtered

python ../extract_sequences-1.py lincRNA_non_redundant_filtered.genes.only_filtered lincRNA_non_redundant_filtered.genes.fa lincRNA_non_redundant_filtered.genes.only_filtered.fa

python ../fasta_header_rename.py lincRNA_non_redundant_filtered.genes.only_filtered.fa lincRNA_final_transcripts.fa
