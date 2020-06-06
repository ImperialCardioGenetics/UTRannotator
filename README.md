# UTRannotator  
A VEP Plugin to annotate high-impact five prime UTR variants either creating new upstream ORFs or disrupting existing upstream ORFs  
  
Currently, it will annotate whether a small variation (1-5bp) including SNVs, indels and MNVs in 5'UTR would have any of the following molecular consequences:  
  
 - [uAUG_gained](#uaug-gained): creating a new start codon AUG  
 - [uAUG_lost](#uaug-lost): removing an existing start codon AUG  
 - [uSTOP_lost](#ustop-lost): removing the stop codon of an existing upstream ORF  
 - [uFrameShift](#uframeshift): creating a frameshift mutation in an existing upstream ORF   
 

# Content 
 [Citation](#citation)
 
 [Requirements](#requirements)
 
 [Installation](#installation)
 
 [Usage](#usage)
 
 [Annotation output](#annotation-output)
 
 [Caveats](#caveats)
 
# Citation  
  
About the role of 5'UTR variants in human genetic disease: 

Whiffin, N., Karczewski, K.J., Zhang, X. et al. Characterising the loss-of-function impact of 5’ untranslated region variants in 15,708 individuals. Nat Commun 11, 2523 (2020). https://doi.org/10.1038/s41467-019-10717-9  

To use our tool, please cite:

Annotating high-impact 5'untranslated region variants with the UTRannotator
Zhang, X., Wakeling, M.N., Ware, J.S, Whiffin,N. bioRxiv 2020.06.03.132266; doi: https://doi.org/10.1101/2020.06.03.132266
    
# Requirements  
- VEP (tested on release-99/202001 and release-100/202005)  
- PERL (tested on version 5.26.2)  
  
# Installation  
To use the plugin with VEP, you would need to add the plugin module in Perl's library path. To do this, you could either:   
  
(1) download all the files of this repository to the VEP default path `$HOME/.vep/Plugins` or  
  
(2) download the repository and add its path to environment variable `$PERL5LIB`.   
  
e.g. Add this line `export PERL5LIB=$PERL5LIB:/path/to/UTRannotator` to `~/.bash_profile`.  
  
# Usage  
  
## Basic Usage  
  
To run the plugin with VEP, you could the following command line:    
  
`vep -i test.vcf --tab -plugin UTRannotator -o test.output`  
  
Currently, the output format is designed for tab-delimited output. 

To be noticed, it's necessary to add option `--minimal` to transform the alleles into minimal representations if it hasn't been transformed beforehand. We have found that this option is necessary especially for variants represented with rs IDs from dbSNP.   
  
## Optional Usage  
  
The plugin could also check whether an input variant disrupts a verified translated uORF.  
  
To use this option, users would pass an evidence file of a list of verified translated uORFs as input.   
  
For translated small ORFs in human, we have curated a list of uORFs previously identified with ribosome profiling from the online repository of small ORFs (www.sorfs.org)  
  
This list is available in the repository:   
  
Genome build GRCh37: `uORF_starts_ends_GRCh37_PUBLIC.txt`  
  
Genome build GRCh38: `uORF_starts_ends_GRCh38_PUBLIC.txt`  

The command to use the file is 

`vep -i test.vcf --tab -plugin UTRannotator,/path/to/uORF_starts_ends_GRCh37_PUBLIC.txt -o test.output`
  
To use a customized list of translated uORF, users would curate a tab-delimited txt file with the following columns:  
  
For example:  
  
`CHR    START_POS GENE    STRAND  TYPE    STOP_POS`  
  
`19  45971469    FOSB    forward five_prime_utr  45971714`  

`START_POS` and `STOP_POS` are the start genomic position and end genomics position of a small ORF respectively. 

  
# Annotation Output  
  
The output annotation from the plugin includes 5 columns:   
  
For any 5'UTR variants, the plugin will first output the number of existing subtype uORFs in the 5'UTR  
  
**Column 1** - existing_InFrame_oORFs: The number of existing inframe overlapping ORFs (inFrame_oORF) already within the 5 prime UTR 

**Column 2** - existing_OutOfFrame_oORFs: The number of existing out-of-frame overlapping ORFs (OutOfFrame_oORF) already within the 5 prime UTR 

**Column 3** - existing_uORFs: The number of existing uORFs with a stop codon within the 5 prime UTR  `

If this 5'UTR is uORF-perturbing, the plugin will output the consequence and detailed annotation of each consequence:   
  
**Column 4** - five_prime_UTR_variant_consequence: Output the variant consequences of a given 5 prime UTR variant: uAUG_gained, uAUG_lost, uSTOP_lost, uFrameshift

**Column 5** - five_prime_UTR_variant_annotation: Output the annotation of a given 5 prime UTR variant 

If a 5'UTR variant perturbs multiple uORFs, the output for each uORF will be concatenated with a vertical bar `|`;   
  
## The detailed annotation for each consequence  
  
 ### uAUG gained
| Annotations                 | Data type | Description                                                                                                                                                                                                        |
|-----------------------------|-----------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| uAUG_gained_type            | String    | The type of of 5’ UTR ORF created, described by one of the following: uORF(with a stop codon in 5’UTR), inframe_oORF (inframe and overlapping  with CDS),OutOfFrame_oORF (out of frame and overlapping with CDS)   |
| uAUG_gained_KozakContext    | String    | The Kozak context sequence of the gained uAUG                                                                                                                                                                      |
| uAUG_gained_KozakStrength   | String    | The Kozak strength of the gained uAUG, described by one of the following values: Weak, Medium or Strong.                                                                                                           |
| uAUG_gained_DistanceToCDS   | Integer   | The distance (number of nucleotides) between the gained  uAUG to CDS                                                                                                                                               |
| uAUG_gained_DistanceFromCap | Integer   | The distance (number of nucleotides) between the gained uAUG to the start of 5’UTR                                                                                                                                 |
| uAUG_gained_DistanceToSTOP  | Integer   | The distance (number of nucleotides) between the gained uAUG to STOP codon (scanning through both the 5’UTR and its downstream CDS). If there is no STOP codon found, it would output NA.                          |

### uAUG lost

| Annotations              | Data type | Description                                                                                                                                     |
|--------------------------|-----------|-------------------------------------------------------------------------------------------------------------------------------------------------|
| uAUG_lost_type           | String    | The type of 5’ UTR ORF lost, described by one of the following: uORF,  inframe_oORF or OutOfFrame_oORF                                          |
| uAUG_lost_KozakContext   | String    | The Kozak context sequence of the lost uAUG                                                                                                     |
| uAUG_lost_KozakStrength  | String    | The Kozak strength of the lost uAUG, described by one of the following values: Weak, Medium or Strong.                                          |
| uAUG_lost_DistanceToCDS  | Integer   | The distance (number of nucleotides) between the lost uAUG to CDS                                                                               |
| uAUG_lost_DistanceToSTOP | Integer   | The distance (number of nucleotides) between the lost uAUG to the nearest stop codon (scanning through both the 5’UTR and its downstream CDS).  |
| uAUG_lost_evidence       | Boolean   | Whether the uORF disrupted by the lost uAUG has any translation evidence. Output NA if no evidence file provided                                                                      |

### uSTOP lost

| Annotations                     | Data type | Description                                                                                                   |
|---------------------------------|-----------|---------------------------------------------------------------------------------------------------------------|
| uSTOP_lost_AltStop              | String    | Whether there is an alternative stop codon downstream within 5’ UTR                                           |
| uSTOP_lost_AltStopDistanceToCDS | Integer   | The distance between the alternative stop codon (if exists) and CDS                                           |
| uSTOP_lost_KozakContext         | String    | The Kozak context sequence of the disrupted uORF                                                              |
| uSTOP_lost_KozakStrength        | String    | The Kozak strength of the disrupted uORF, described by one of the following values: Weak, Medium or Strong.   |
| uSTOP_lost_FrameWithCDS         | String    | The frame of the uORF with respect to CDS, described by inFrame or outOfFrame.                                |
| uSTOP_lost_evidence             | Boolean   | Whether the uORF disrupted by the lost stop codon has any translation evidence. Output NA if no evidence file provided.                               |

### uFrameShift

| Annotations                     | Data type | Description                                                                                                        |
|---------------------------------|-----------|--------------------------------------------------------------------------------------------------------------------|
| uFrameshift_ref_type            | String    | The type of uORF with the reference allele, described by one of following: uORF, inframe_oORF or OutOfFrame_oORF   |
| uFrameshift_StartDistanceToCDS  | Integer   | The distance between the start codon of the disrupting uORF and CDS                                                |
| uFrameshift_alt_type            | String    | The type of uORF with the alternative allele, described by one of following: uORF, inframe_oORF or OutOfFrame_oORF |
| uFrameshift_KozakContext        | String    | The Kozak context sequence of the disrupted uORF                                                                   |
| uFrameshift_KozakStrength       | String    | The Kozak strength of the disrupted uORF, described by one of the following values: Weak, Medium or Strong.        |
| uFrameshift_evidence            | Boolean   | Whether the disrupted uORF has any translation evidence. Output NA if no evidence file provided                                                          |


# Caveats 
- only tested on SNVs, small indels (1-5bp) and MNV (1-5bp)  
- only consider canonical start codon  
