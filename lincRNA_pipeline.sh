#!/bin/bash
while getopts ":b:c:hr:" opt; do
  case $opt in
    b)
      blastfile=$OPTARG
      ;;
    c)
      comparefile=$OPTARG
      ;;
    h)
      echo "USAGE : test.sh -c cuffcompare_output -r reference_genome -b blast_file"
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
