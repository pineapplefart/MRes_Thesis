makeblastdb -in assembly.fasta -dbtype nucl
blastn \
  -perc_identity 95 \
  -qcov_hsp_perc 90 \
  -query promoter_set_cds_feb20.fasta \
  -db assembly.fasta \
  -out cds_vs_genome.tsv \
  -outfmt "6 qseqid sseqid sstart send sstrand pident length" \
  -max_target_seqs 1 \
  -evalue 1e-20