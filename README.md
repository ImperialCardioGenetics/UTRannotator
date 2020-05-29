# 5'UTR annotator
VEP Plugin to annotate high-impact five prime UTR variants either creating new upstream ORFs or disrupting existing upstream ORFs

Currently, it will annotate whether a small variation (1-5bp) in 5'UTR would have any of the following molecular consequences:

 - uAUG-gained: creating a new start codon AUG
 - uAUG-lost: removing an existing start codon AUG
 - uSTOP-lost: removing the stop codon of an existing upstream ORF
 - uFrameShift: creating a frameshift mutation in an existing upstream ORF 

# Citation

Whiffin, N., Karczewski, K.J., Zhang, X. et al. Characterising the loss-of-function impact of 5’ untranslated region variants in 15,708 individuals. Nat Commun 11, 2523 (2020). https://doi.org/10.1038/s41467-019-10717-9


# Caveats 
- only tested on SNVs, small indels (1-5bp) and MNV (1-5bp)
- only consider canonical start codon

# Requirements
- VEP (tested on release-99/202001)
- PERL (tested on version 5.26.2)

# Usage 
To use the plugin with VEP, you would need to add the plugin module in Perl's library path. To do this, you could either: 

(1) copy all the files of this repository to the VEP default path `$HOME/.vep/Plugins` or
(2) copy the repository and add its path to environment variable `$PERL5LIB`. 

e.g. Add this line `export PERL5LIB=$PERL5LIB:/path/to/5primeUTRannotator` to `~/.bash_profile`.

The Plugin could be run with VEP using the following command (if using hg19 genome build): 

`vep -i test.vcf --database --hgvs --tab --port 3337 --minimal -plugin five_prime_UTR_annotator,/path/to/uORF_starts_ends_GRCh37_PUBLIC.txt -o test.output`. 

To be noticed, it's necessary to add option `--minimal` to transform the alleles into minimal representations if it hasn't been transformed beforehand. We have found that this option is necessary especially for variants represented with rs IDs from dbSNP. 


