use experimental 'smartmatch';

sub uFrameshift{
    #assumption: there is no insertion or deletion happening before stop codon. If there are, will output as frameshift or inframe indels.


    #Description: annotate if a five_prime_UTR_varint removes a stop codon of an existing uORF(given that uORF doesn't not change)


    #Returntype: hashref

    my ($self, $variant_info, $UTR_info) = @_;
    my %flip;
    $flip{'A'}='T';
    $flip{'C'}='G';
    $flip{'G'}='C';
    $flip{'T'}='A';

    my %kozak_strength;
    $kozak_strength{1}='Weak';
    $kozak_strength{2}='Moderate';
    $kozak_strength{3}='Strong';

    my $chr = $variant_info->{chr};
    my $pos = $variant_info->{pos};
    my $ref = $variant_info->{ref};
    my $alt = $variant_info->{alt};

    my @sequence = split //, $UTR_info->{seq};
    my $strand = $UTR_info->{strand};
    my $utr_length = @sequence;

    my %existing_uORF = {};
    %existing_uORF = %{$self->existing_uORF(\@sequence)};

    #return annotators
	my $uFrameshift_ref_type = ""; # the type of uORF with the reference allele
    my $uFrameshift_StartDistanceToCDS = ""; # the distance between the start codon of the disrupted uORF and CDS
    my $uFrameshift_alt_type = ""; # the type of uORF with the alternative allele
    my $uFrameshift_KozakContext = ""; # the Kozak context sequence of the disrupted uORF
    my $uFrameshift_KozakStrength = ""; # the Kozak strength of the disrupted uORF
    my $uFrameshift_evidence = ""; # whehter there is translation evidence for the disrupted uORF

    #indicate whether the variant ever introduce a frameshift variant
    my $flag_uORF;
    my $output_flag = "";

    my $current_kozak;
    my $current_kozak_strength;

    my %result = (); #result is a hash table with two elements: $flag and $output_effects
    my $output_effects="";

    my $ref_coding = $self->get_ref_coding($ref);
    my $alt_coding = $self->get_alt_coding($alt,$strand);

	#skip alleles with same length
	return {} unless(length($ref_coding) ne length($alt_coding));

 	my ($mut_pos, $end_pos) = $self->get_allele_exon_pos($strand, $pos, $ref_coding, $UTR_info);
    return {} unless(defined($mut_pos)&defined($end_pos));

    #if it's a deletion at the boundary of exon and intron, we would skip the annotation

	my $mut_utr_seq = $self->mut_utr_sequence(\@sequence,$mut_pos,$ref_coding,$alt_coding,$strand);
  	my @mut_utr_seq = split //,$mut_utr_seq;
    my $length = @mut_utr_seq;

    my @start = sort {$a <=> $b} (keys %existing_uORF);
    my $offset;


	#check for each uORF
    if(%existing_uORF){
    for (my $i=0;$i<@start;$i++){
        $flag_uORF=0;
        my $start_pos = $start[$i];
        my @stops = sort {$a <=> $b} @{$existing_uORF{$start_pos}};
        my $stop_pos=$stops[0];


		#if it's entirely in the uORF
		if(($mut_pos>=$start_pos)&($end_pos<=$stop_pos+2)){
		#snp and insertion could only be annotated here
			$flag_uORF = abs(length($ref_coding)-length($alt_coding))%3;
		}elsif(($mut_pos<$start_pos)&($end_pos>=$start_pos)){

			$offset = $start_pos-$mut_pos;
			#if it's an deletion
			if(length($alt_coding) eq 0){
			#Given that (1)length($ref_coding) = $end_pos-$start_pos
			#(2)length($ref_coding)> = $offset
			$flag_uORF = abs(length($ref_coding)-$offset)%3;}
			elsif(length($alt_coding)>=$offset){
			$flag_uORF = abs((length($ref_coding)-$offset)-(length($alt_coding)-($offset)))%3;
			}

		}

        if($flag_uORF){

                #getting the Kozak context and Kozak strength of the start codon

  				if ((($start_pos-3)>=0)&&($sequence[($start_pos+3)])){
  					$current_kozak = $sequence[($start_pos-3)].$sequence[($start_pos-2)].$sequence[$start_pos-1]."ATG".$sequence[$start_pos+3];
  				}
  				else{
  					$current_kozak = '-';
  				}

                if ($current_kozak !~ /-/){
                    my @split_kozak = split //, $current_kozak;
                    $current_kozak_strength = 1;
                    if ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))&&($split_kozak[6] eq 'G')){
                        $current_kozak_strength = 3;
                    }
                    elsif ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))||($split_kozak[6] eq 'G')){
                        $current_kozak_strength = 2;
                    }
                }

                $uFrameshift_KozakContext=$current_kozak;
                $uFrameshift_KozakStrength=$kozak_strength{$current_kozak_strength}? $kozak_strength{$current_kozak_strength}:$current_kozak_strength;


                #the annotation of the original uORF
                if (@{$existing_uORF{$start_pos}}){ #if there is stop codon within 5'UTR
                    $uFrameshift_ref_type = "uORF"
                }elsif (($utr_length-$start_pos) % 3){
                	$uFrameshift_ref_type = "OutOfFrame_oORF";
                    }
                else{
                    $uFrameshift_ref_type = "inFrame_oORF";
                    }

                $uFrameshift_StartDistanceToCDS = $utr_length - $start_pos;

                #if there is an alternative stop codon in the mutant uORF sequence
                my %mut_uORF = %{$self->existing_uORF(\@mut_utr_seq)};

                # the sequence before mut_pos should not be changed. Thus start_pos shall still correspond to a start codon;
                # if the sequence is indeed very short as such ATGTGA

                my @mut_stops;
				#only if the start codon hasn't been affected
                if(exists($mut_uORF{$start_pos})){
                @mut_stops = sort {$a <=> $b} @{$mut_uORF{$start_pos}}
                };

                if (@mut_stops>0){
                $uFrameshift_alt_type = "uORF";
                } #if there is no alternative stop codon
                elsif(($length-$start_pos)%3) {
                    $uFrameshift_alt_type = "OutOfFrame_oORF";
                }
                else{
                    $uFrameshift_alt_type = "inFrame_oORF";
                }

            	#find evidence from sorf
                my %utr_pos = %{$self->chr_position($UTR_info)};
                my $start_chr_pos = $utr_pos{$start_pos};
            	#find evidence from sorf
            	##TODO: fix the finding evidence of uORF
            	my $query = ($chr=~/chr/i)?$chr.":".$start_chr_pos:"chr".$chr.":".$start_chr_pos;
            	$uFrameshift_evidence=$self->{uORF_evidence}->{$query}?"True":"False";


                my %uORF_effect = (
                "uFrameShift_ref_type" => $uFrameshift_ref_type,
                "uFrameShift_ref_uAUGDistanceToCDS" => $uFrameshift_StartDistanceToCDS,
                "uFrameShift_alt_type" => $uFrameshift_alt_type,
                "uFrameShift_KozakContext" => $uFrameshift_KozakContext,
                "uFrameShift_KozakStrength" => $uFrameshift_KozakStrength,
                "uFrameShift_Evidence" => $uFrameshift_evidence,
	             );

                $output_flag = $output_flag? $output_flag."|"."uFrameShift":"uFrameShift";
                $output_effects = $output_effects? $output_effects."|".$self->transform_hash_to_string(\%uORF_effect):$self->transform_hash_to_string(\%uORF_effect);
            }

    }
    }

    $result{'uFrameShift_flag'} = $output_flag;
    $result{'uFrameShift_effect'} = $output_effects;

    return \%result;
}

1;