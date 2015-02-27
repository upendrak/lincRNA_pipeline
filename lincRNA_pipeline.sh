#!/bin/bash
if [ "$1" == "--help" ]
then
    echo "Usage: sh $0 cuffcompare_output reference_genome blast_file"
    exit
else
wget -nv https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz && tar xvf 2.0.1 && rm -r 2.0.1
makeblastdb -in $3 -dbtype nucl -out $3.blast.out
grep '"u"' $1 | \
gffread -w transcripts_u.fa -g $2 - && \
python2.7 get_gene_length_filter.py transcripts_u.fa transcripts_u_filter.fa && \
TransDecoder-2.0.1/TransDecoder.LongOrfs -t transcripts_u_filter.fa
sed 's/ .*//' transcripts_u_filter.fa | grep ">" | sed 's/>//' > transcripts_u_filter.fa.genes
cd transcripts_u_filter.fa.transdecoder_dir
sed 's/|.*//' longest_orfs.cds | grep ">" | sed 's/>//' | uniq > longest_orfs.cds.genes
grep -v -f longest_orfs.cds.genes ../transcripts_u_filter.fa.genes > longest_orfs.cds.genes.not.genes
sed 's/^/>/' longest_orfs.cds.genes.not.genes > temp && mv temp longest_orfs.cds.genes.not.genes
python ../extract_sequences.py longest_orfs.cds.genes.not.genes ../transcripts_u_filter.fa longest_orfs.cds.genes.not.genes.fa
blastn -query longest_orfs.cds.genes.not.genes.fa -db ../$3.blast.out -out longest_orfs.cds.genes.not.genes.fa.blast.out -outfmt 6
python ../filter_sequences.py longest_orfs.cds.genes.not.genes.fa.blast.out longest_orfs.cds.genes.not.genes.fa.blast.out.filtered
grep -v -f longest_orfs.cds.genes.not.genes.fa.blast.out.filtered longest_orfs.cds.genes.not.genes.fa > lincRNA_final.fa
fi