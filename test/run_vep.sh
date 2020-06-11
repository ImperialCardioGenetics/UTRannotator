vep -i test.vcf --database --port 3337 --plugin UTRannotator,../uORF_starts_ends_GRCh37_PUBLIC.txt --tab --force_overwrite -o test.output
