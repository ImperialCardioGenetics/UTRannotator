vep -i test_grch37.vcf --database --species homo_sapiens --port 3337 --plugin UTRannotator,../uORF_starts_ends_GRCh37_PUBLIC.txt --tab --force_overwrite -o test_output_grch37.txt
vep -i test_grch38.vcf --database --species homo_sapiens --assembly GRCh38 --plugin UTRannotator,../uORF_starts_ends_GRCh38_PUBLIC.txt --tab --force_overwrite -o test_output_grch38.txt
vep -i test_mus.vcf --database --species mus_musculus --assembly GRCm38 --plugin UTRannotator --tab --force_overwrite -o test_mus.output
