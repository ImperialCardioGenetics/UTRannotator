# uORF_annotator
VEP Plugin to annotate high-impact five prime UTR variants either creating new upstream ORFs or disrupting existing upstream ORFs

# Caveats 
- only tested on SNVs, small indels (1-5bp) and MNV (1-5bp)
- only consider canonical start codon

# Usage 
To use the plugin with VEP, you would need to add the plugin module in Perl's library path. To do this, you could either 
(1) copy all the files of this repository to the VEP default path `~/.vep/Plugins` or
(2) copy the repository and add its path to environment variable `$PERL5LIB`.

The Plugin could be run with VEP using the following command (if using hg19 genome build): 

`vep -i test_var.vcf --database --hgvs --tab --port 3337 --force_overwrite -plugin five_prime_UTR_annotator,/path/to/uORF_starts_ends_GRCh37_PUBLIC.txt -o test_var.output`
