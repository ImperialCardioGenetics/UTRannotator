# uORF_annotator
VEP Plugin to annotate high-impact five prime UTR variants

# Caveats 
- only tested on SNVs, small indels (1-5bp) and MNV (1-5bp)
- only consider canonical start codon

# Usage 
The plugin is applied with VEP using the following command (if using hg19 genome build): 

`vep -i test_var.vcf --database --hgvs --tab --port 3337 --force_overwrite -plugin five_prime_UTR_annotator,./uORF_starts_ends_GRCh37_PUBLIC.txt -o test_var.output
