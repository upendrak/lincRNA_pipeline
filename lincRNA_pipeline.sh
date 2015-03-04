#!/bin/bash
while getopts ":b:c:hr:" opt; do
  case $opt in
    b)
      blastfile=$OPTARG
      ;;
    c)
      comparefile=$OPTARG
      ;;
    h)	echo 	"USAGE : sh lincRNA_pipeline.sh 
		          -c 	</path/to/cuffcompare_output file>
	 		  -r 	</path/to/reference_genome file>
			  -b 	</path/to/blast_file>"
      exit 1
      ;;
    r)
      referencegenome=$OPTARG
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

blastn -query longest_orfs.cds.genes.not.genes.fa -db ../$blastfile.blast.out -out longest_orfs.cds.genes.not.genes.fa.blast.out -outfmt 6

python ../filter_sequences.py longest_orfs.cds.genes.not.genes.fa.blast.out longest_orfs.cds.genes.not.genes.fa.blast.out.filtered

grep -v -f longest_orfs.cds.genes.not.genes.fa.blast.out.filtered longest_orfs.cds.genes.not.genes.fa \
	> lincRNA_final.fa

cap3 lincRNA_final.fa

cat lincRNA_final.fa.cap.singlets lincRNA_final.fa.cap.contigs \
  > lincRNA_final_merged.fa

python ../fasta_header_rename.py lincRNA_final_merged.fa lincRNA_final_merged_renamed.fa
