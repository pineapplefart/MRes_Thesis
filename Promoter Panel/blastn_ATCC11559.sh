makeblastdb -in GCA_046865755.1/cds_from_genomic.fna -dbtype nucl -out target_strain_db
blastn -query promoter_set_cds_feb20.fasta -db target_strain_db -out results.out -outfmt 6