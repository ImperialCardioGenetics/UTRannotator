

sub uSTOP_gained{
	#Description: annotate if a five_prime_UTR_variant creates ATG
	#Returntype: hashref

	my ($self, $variant_info,$UTR_info) = @_;
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

	#return annotators
	my $uSTOP_gained_DistanceToCDS = "";  # the distance between the gained uAUG to CDS
	my $uSTOP_gained_KozakContext = "";  # the Kozak context sequence of the gained uAUG
	my $uSTOP_gained_KozakStrength = ""; # the Kozak strength of the gained uAUG
	my $uSTOP_gained_ref_type = ""; # the type of uORF created - any of the following: uORF, inframe_oORF,OutOfFrame_oORF
	my $uSTOP_gained_DistanceFromCap = ""; # the distance between the gained uAUG to the start of the five prime UTR
#    my $uSTOP_gained_DistanceToStop = ""; #the distance between the gained uAUG to stop codon (could be in CDS)

	#indicate whether the variant creates a ATG
	my $flag = 0;
	my $current_kozak = "";
	my $current_kozak_strength ="";

	#the relative position of input variant in the UTR sequence

    my %result = (); #result is a hash table with two elements: $flag and $output_effects
    my $output_flag = 0;
    my $output_effects = "";

    my $ref_coding = $self->get_ref_coding($ref);
    my $alt_coding = $self->get_alt_coding($alt,$strand);

	my ($mut_pos, $end_pos) = $self->get_allele_exon_pos($strand, $pos, $ref_coding, $UTR_info);
    return {} unless(defined($mut_pos)&defined($end_pos));

	my $mut_utr_seq = $self->mut_utr_sequence(\@sequence,$mut_pos,$ref_coding,$alt_coding,$strand);
  	my @mut_utr_seq = split //,$mut_utr_seq;
    my $mut_utr_length = @mut_utr_seq;

    my @start = @{$self->get_ATG_pos(\@sequence)};


    if(@start){
    for (my $i=0;$i<@start;$i++){
        $flag_uORF=0;
        my $start_pos = $start[$i];

        # to set the range for searching new STOP codon: not including start_codon, (1) for uORFs - start_codon...stop_codon (2) for oORFs - start_codon...end of 5'UTR
        my $check_point = $start_pos;
        # the checking end point of eligible area
        #For uORF: start_pos .. check_point
        #For overlapping ORF: start_pos .. 3' end of 5'UTR sequence

        #if the variant is entirely in the uORF (within 5'UTR
        if(exists($existing_uORF{$start_pos})) {
            my @stops = sort {$a <=> $b} @{$existing_uORF{$start_pos}};
            $check_point = $stops[0]-1;;
        }else{
        #if the existing uORF is an oORF
            $check_point = $utr_length-1;
        }

        # only check the ones between start and stop codons and
        # ignore the cases at the boundary of start and stop codon since the effect on start/stop would be evaluated anyway
		if(($mut_pos>=$start_pos+3)&($end_pos<=$check_point)){
		#snp and insertion could only be annotated here

        #for snps and deletion
        if (length($ref_coding)){

        if ($mut_pos+length($ref_coding)-1<$stop_pos) {next;}

        } elsif ($mut_pos<$stop_pos) {next;}   #for insertion


        #for deletion, it definitely disrupting the stop codon.

        if (length($alt_coding) eq 0) {
        $flag_uORF=1;
        }else{
        my $mut_codon = @mut_utr_seq[$stop_pos].@mut_utr_seq[$stop_pos+1].@mut_utr_seq[$stop_pos+2];
        if ($mut_codon ~~ @stop_codons) {next;}
        $flag_uORF=1;
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
                    $uFrameshift_ref_type = "InFrame_oORF";
                    }

                $uFrameshift_StartDistanceToCDS = $utr_length - $start_pos;

                #if there is an alternative stop codon in the mutant uORF sequence
                my %mut_uORF = %{$self->existing_uORF(\@mut_utr_seq)};

                my @mut_stops;
				#Since we only evaluate frameshift within uORFs, thus the start codon shall not be affected.
                #Thus there shall always be an uORF existed after mutation.
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
                    $uFrameshift_alt_type = "InFrame_oORF";
                }

            	#find evidence from sorfs.org
                my %utr_pos = %{$self->chr_position($UTR_info)};
                my $start_chr_pos = $utr_pos{$start_pos};
            	#find evidence from sorf
            	##TODO: fix the finding evidence of uORF
            	my $query = ($chr=~/chr/i)?$chr.":".$start_chr_pos:"chr".$chr.":".$start_chr_pos;
                $uFrameshift_evidence=$self->{uORF_evidence}->{$query}?"True":"False";

                my %uORF_effect = (
                "uFrameShift_ref_type" => $uFrameshift_ref_type,
                "uFrameShift_ref_StartDistanceToCDS" => $uFrameshift_StartDistanceToCDS,
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




    my @mut_seq = @mut_utr_seq[$mut_pos-2..$mut_pos+length($alt_coding)+1];
    my @mut_atg_pos = @{$self->get_ATG_pos(\@mut_seq)};

    # check whether there is a ATG codon is to check the length of the array
    #if there is a ATG codon, flag it and output the relative positive of A
    my $pos_A;
    if (@mut_atg_pos){
    #get the pos of A (in ATG) in the UTR sequence
    $pos_A=$mut_pos-2+$mut_atg_pos[0];
    $flag=1;
    }


    if ($flag ==1){

  ################################################################################
  #annotator 1: get the distance to the start codon of the main ORF
  ################################################################################
  		$uAUG_gained_DistanceToCDS = $mut_utr_length-$pos_A;
  		$uAUG_gained_DistanceFromCap = $pos_A;
  ################################################################################
  #annotator 2: determine kozak context;
  ################################################################################
  		if ((($pos_A-3)>=0)&&($mut_utr_seq[($pos_A+3)])){
  			$current_kozak = $mut_utr_seq[($pos_A-3)].$mut_utr_seq[($pos_A-2)].$mut_utr_seq[$pos_A-1]."ATG".$mut_utr_seq[$pos_A+3];
  		}
  		else{
  			$current_kozak = '-';
  		}
  		#get the strength of kozak context
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



		$uAUG_gained_KozakContext=$current_kozak;
        $uAUG_gained_KozakStrength=$kozak_strength{$current_kozak_strength}? $kozak_strength{$current_kozak_strength}:$current_kozak_strength;


  ################################################################################
  #annotator 3: Type of new ORF with respect the main ORF: uORF/Overlapping_inFrame/Overlapping_outofframe
  ################################################################################

        #check what kind of uORF does that correspond to?
        #first check whether it's overlapping with CDS
        my %existing_utr_uorf = %{$self->existing_uORF(\@mut_utr_seq)};


        if(exists($existing_utr_uorf{$pos_A})){ #if there is stop codon within 5'UTR
             $uAUG_gained_type = "uORF";
        }else{

            if(($mut_utr_length - $pos_A) % 3){
                $uAUG_gained_type = "OutOfFrame_oORF";
            }else{
                $uAUG_gained_type = "inFrame_oORF";
            }
        }

        my @overlapping_seq = split //, $mut_utr_seq.$UTR_info->{cds_seq};
        my %existing_uORF = %{$self->existing_uORF(\@overlapping_seq)};
        if(exists($existing_uORF{$pos_A})){
        my @stop_pos_array = sort{$a<=>$b}@{$existing_uORF{$pos_A}};
        my $stop_pos = @stop_pos_array[0];
        $uAUG_gained_DistanceToStop = $stop_pos-$pos_A;
        }else{
        $uAUG_gained_DistanceToStop = "NA";
        }


  		my %uORF_effect = (
        "uAUG_gained_KozakContext" => $uAUG_gained_KozakContext,
        "uAUG_gained_KozakStrength" => $uAUG_gained_KozakStrength,
        "uAUG_gained_DistanceToCDS" => $uAUG_gained_DistanceToCDS,
        "uAUG_gained_type" => $uAUG_gained_type,
        "uAUG_gained_DistanceToStop" => $uAUG_gained_DistanceToStop,
        "uAUG_gained_DistanceFromCap" => $uAUG_gained_DistanceFromCap,
      );

       $output_flag = "uAUG_gained";
       $output_effects = $self->transform_hash_to_string(\%uORF_effect);


  	}
    $result{'uAUG_gained_flag'} = $output_flag;
    $result{'uAUG_gained_effect'} = $output_effects;

    return \%result;

}

1;
