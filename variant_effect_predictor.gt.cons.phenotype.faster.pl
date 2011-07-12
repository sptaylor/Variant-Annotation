#!/usr/bin/perl


##NEW June 1 2011, Added Allele Frequencies
=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Variant Effect Predictor - a script to predict the consequences of genomic variants

Version 2.0

by Will McLaren (wm2@ebi.ac.uk)
=cut

use strict;
use Getopt::Long;
use FileHandle;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::TranscriptVariationAdaptor;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(MAX_DISTANCE_FROM_TRANSCRIPT);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(unambiguity_code);
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::Seq;
use Bio::Matrix::IO;
use Bio::AlignIO;
use Bio::DB::EUtilities;
use Bio::SeqIO;
use SWISS::Entry;	    	
use SWISS::CCs;
use SWISS::DEs;


	
# configure from command line opts
my $config = &configure(scalar @ARGV);

# run the main sub routine
&main($config);

# this is the main sub-routine - it needs the configured $config hash
sub main {
	my $config = shift;
	
	debug("Starting...") if defined $config->{verbose};
	
	my $species = $config->{species};

	# get adaptors
	$config->{vfa} = $config->{reg}->get_adaptor($species, 'variation', 'variationfeature');
	$config->{va} = $config->{reg}->get_adaptor($species, 'variation', 'variation');
	$config->{vaa} = $config->{reg}->get_adaptor($species, 'variation', 'variationannotation');
	$config->{tva} = $config->{reg}->get_adaptor($species, 'variation', 'transcriptvariation');
	$config->{ma} = $config->{reg}->get_adaptor( 'Multi', 'compara', 'Member' );
	$config->{ha} = $config->{reg}->get_adaptor( 'Multi', 'compara', 'Homology' );
	#$config->{fa} = $config->{reg}->get_adaptor( 'Multi', 'compara', 'Family' );
	$config->{goa} = $config->{reg}->get_adaptor( 'Multi', 'Ontology', 'GOTerm' );
	$config->{sa} = $config->{reg}->get_adaptor( $species, 'core', 'slice' );
	$config->{ta} = $config->{reg}->get_adaptor( $species, 'core', 'transcript' );
	$config->{ga} = $config->{reg}->get_adaptor( $species, 'core', 'gene' );
	$config->{translation_adaptor} = $config->{reg}->get_adaptor( $species, 'core', 'translation' );
	#$config->{genomic_align_block_adaptor} = $config->{reg}->get_adaptor( 'Multi', 'compara', 'GenomicAlignBlock');
	$config->{pga} = $config->{reg}->get_adaptor("human","variation","populationgenotype");  
	$config->{pa} = $config->{reg}->get_adaptor("human","variation","population");
	
	#For get_conservation_scores
	$config->{mlss_adaptor} = $config->{reg}->get_adaptor( 'Multi', 'compara', 'MethodLinkSpeciesSet');
	$config->{genome_db_adaptor} = $config->{reg}->get_adaptor('Multi', 'compara', "GenomeDB");
	
	# get fake ones for species with no var DB
	if(!defined($config->{vfa})) {
		$config->{vfa} = Bio::EnsEMBL::Variation::DBSQL::VariationFeatureAdaptor->new_fake($species);
	}
	else {
		$config->{vfa}->db->include_failed_variations($config->{include_failed}) if defined($config->{vfa}->db) && $config->{vfa}->db->can('include_failed_variations');
	}
	
	# check we got slice adaptor - can't continue without a core DB
	die("ERROR: Could not connect to core database\n") unless defined $config->{sa} and defined $config->{ga};
	
	
	# create a hash to hold slices so we don't get the same one twice
	my %slice_hash = ();
	my @new_vfs;
	my @ends;
	my %vf_hash;

	
	my $transcript_cache;
	
	my $line_number = 0;
	my $vf_count;
	my $in_file_handle = $config->{in_file_handle};
	my $time = &get_time;
	print $time;
	
	# read the file
	while(<$in_file_handle>) {
	  chomp;
	  
	  $line_number++;
	  
	  # header line?
	  next if /^\#/;
	  
	  # some lines (pileup) may actually parse out into more than one variant)
	  foreach my $sub_line(@{&parse_line($config, $_)}) {
		
		# get the sub-line into named variables
		my ($chr, $start, $end, $allele_string, $strand, $var_name, $genotypes) = @{$sub_line};
		
		# non-variant line from VCF
		next if $chr eq 'non-variant';
		
		# fix inputs
		$chr =~ s/chr//ig;
		$strand = ($strand =~ /\-/ ? "-1" : "1");
		$allele_string =~ tr/acgt/ACGT/;
		
		# sanity checks
		unless($start =~ /^\d+$/ && $end =~ /^\d+$/) {
		  warn("WARNING: Start $start or end $end coordinate invalid on line $line_number\n") unless defined $config->{quiet};
		  next;
		}
		
		unless($allele_string =~ /([ACGT-]+\/*)+/) {
		  warn("WARNING: Invalid allele string $allele_string on line $line_number\n") unless defined $config->{quiet};
		  next;
		}
		
		# now get the slice
		my $slice;
	   
		# check if we have fetched this slice already
		if(defined $slice_hash{$chr}) {
		  $slice = $slice_hash{$chr};
		}
	   
		# if not create a new one
		else {
		  
		  # first try to get a chromosome
		  eval { $slice = $config->{sa}->fetch_by_region('chromosome', $chr); };
		  
		  # if failed, try to get any seq region
		  if(!defined($slice)) {
			  $slice = $config->{sa}->fetch_by_region(undef, $chr);
		  }
		  
		  # if failed, warn and skip this line
		  if(!defined($slice)) {
			  warn("WARNING: Could not fetch slice named $chr on line $line_number\n") unless defined $config->{quiet};
			  next;
		  }	
		  
		  # store the hash
		  $slice_hash{$chr} = $slice;
		}
		
		# check reference allele if requested
		if(defined $config->{check_ref}) {
		  my $ref_allele = (split /\//, $allele_string)[0];
		  
		  my $ok = 0;
		  my $slice_ref_allele;
		  
		  # insertion, therefore no ref allele to check
		  if($ref_allele eq '-') {
			  $ok = 1;
		  }
		  else {
			  my $slice_ref = $slice->sub_Slice($start, $end, $strand);
			  
			  if(!defined($slice_ref)) {
				  warn "WARNING: Could not fetch sub-slice from $start\-$end\($strand\) on line $line_number" unless defined $config->{quiet};
			  }
			  
			  else {
				  $slice_ref_allele = $slice_ref->seq;
				  $ok = ($slice_ref_allele eq $ref_allele ? 1 : 0);
			  }
		  }
		  
		  if(!$ok) {
			  warn
				  "WARNING: Specified reference allele $ref_allele ",
				  "does not match Ensembl reference allele",
				  ($slice_ref_allele ? " $slice_ref_allele" : ""),
				  " on line $line_number" unless defined $config->{quiet};
			  next;
		  }
		}
	   
		# create a new VariationFeature object
		my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new(
		  -start => $start,
		  -end => $end,
		  -slice => $slice,           # the variation must be attached to a slice
		  -allele_string => $allele_string,
		  -strand => $strand,
		  -map_weight => 1,
		  -adaptor => $config->{vfa},           # we must attach a variation feature adaptor
		  -variation_name => (defined $var_name ? $var_name : $chr.'_'.$start.'_'.$allele_string),
		  -source => $genotypes,
		  
		);
		
		if(defined $config->{whole_genome}) {
			push @{$vf_hash{$chr}{int($start / $config->{chunk_size})}{$start}}, $new_vf;
			#my $what =int($start / $config->{chunk_size});
			$vf_count++;
			
			if($vf_count == $config->{buffer_size}) {
				&print_consequences($config, &whole_genome_fetch($config, \%vf_hash, $transcript_cache), %slice_hash);
				%vf_hash = ();
				$vf_count = 0;
			}
		}
		else {
			&print_consequences($config, [$new_vf]);
			$vf_count++;
			debug("Processed $vf_count variants") if $vf_count =~ /0$/ && defined($config->{verbose});
		}
	  }
	}
	
	# if in whole-genome mode, finish off the rest of the buffer
	if(defined $config->{whole_genome} && %vf_hash) {
		&print_consequences($config, &whole_genome_fetch($config, \%vf_hash, $transcript_cache), %slice_hash);
	}
	
	debug("Finished!") if defined $config->{verbose};
}

sub print_consequences {
	my $config = shift;
	my $vfs = shift;
	my $slice_hash = shift;
	
	
	my ($filter, $qual);
	if ($config->{format}=='vcf') {
			my @data = (split /\s+/, $_);
			$filter=$data[6];
			$qual=$data[5];
	}
	
	my $out_file_handle = $config->{out_file_handle};
	
	my %parsed;
	
	foreach my $new_vf(@$vfs) {
		# find any co-located existing VFs
		my ($existing_vf, $existing_vf_id);
		my $all_freq;
		my @pops;
		my @freqs;
		my @counts;
		my $alt_frequency;
		#$existing_vf = &find_existing($new_vf) if defined $config->{check_existing} && $config->{check_existing} == 1;
		 if (defined $config->{check_existing} && $config->{check_existing} == 1) {
			$existing_vf = &find_existing($new_vf);
			
			unless (!defined ($existing_vf) ) {
			my $v = $config->{va}->fetch_by_name($existing_vf);
			$existing_vf_id=$v->{dbID};
			if (defined ($existing_vf_id) ){
				#Allele frequencies:
			if(defined($new_vf->adaptor->db)) {
				my $sth = $new_vf->adaptor->db->dbc->prepare(qq{
					SELECT distinct 
				      A.allele,
				      A.count,
				      i.population_sample_id,
				      S.name
				      from 
				      allele as A,
				      individual_population as i,
					sample as S
				     where S.sample_id = i.population_sample_id and
				        A.sample_id=i.population_sample_id and
				      A.variation_id= ?
				      });
	
				my ($ref_allele, $alt_allele)=split /\//, $new_vf->allele_string;
				$sth->execute($existing_vf_id);		
						
				my ($allele, $count, $pop_id, $sample_name);
				$sth->bind_columns( \$allele, \$count, \$pop_id, \$sample_name);
		
				#my %by_source;

				#my @frequencies;
				my $ref_idx = 0;
				my $alt_idx = 0;
				while ($sth->fetch) {
					
					if ($allele eq $alt_allele ) {
						$alt_idx+=$count;
					}
					else {
						$ref_idx+=$count;
					}
				}
				my $total = $ref_idx+$alt_idx;
				$alt_frequency=($alt_idx/$total) unless ($total eq 0);
				
			}
			}
				
				my @alleles = @{$v->get_all_Alleles()};
			
				#my @population_genotypes=@{$v->get_all_PopulationGenotypes};
				#foreach my $pg(@population_genotypes){
				#	my $count= $pg->count;
				#	my $frequency=$pg->frequency;
				#	my $pop=$pg->population;
				#}
				my $var_allele= (split /\//, $new_vf->allele_string)[0];
				foreach my $allele(@alleles) {
					if ( ($allele->allele eq $var_allele) && defined($allele->frequency)  && defined($allele->population) ) {
						
						unless ($allele->is_failed() ) {
						$all_freq->{(split /\:/, $allele->population->name)[1]}=($allele->frequency.",".$allele->count());
						my $pop_size=$allele->population->size;
						push @pops, (split /\:/, $allele->population->name)[1];
						push @freqs, $allele->frequency;
						push @counts, $allele->count();
						}
					}
				}
			}
		 }
		

		$existing_vf ||= '-';

		unless (scalar (@{$new_vf->get_all_TranscriptVariations}) ) {
			warn("No Transcript found for variant ".$new_vf->variation_name."\n");
			print $out_file_handle
				  $new_vf->variation_name, "\t",
				  $new_vf->seq_region_name, ":",
				  &format_coords($new_vf->start, $new_vf->end), "\t",
				  (split /\//, $new_vf->allele_string)[0], "\t",
				  (split /\//, $new_vf->allele_string)[1], "\t",
				  $new_vf->source, "\t",
				   ".","\t",
				  ".","\t",
				  ".", "\t",			 
				  ".","\t",
				  ".","\t",
				  "No Transcript \t",
				  ("."), "\t",
				  ".", "\t",
				  $existing_vf, "\t",
				  ".","\t",
				  ".","\t",
				  ".","\t",
				  ".","\t",
				  ".","\t",
				  ".","\t",
				  ".","\t",
				  ".","\t",
				  ".","\t",
				  ".","\t",
				  ".","\t",
				  ".","\t",
				  ".","\t",
				  ".","\t",
				  ".","\n";
		}	  
	
		next if(defined $config->{damaging} && ( int($new_vf->most_severe_OverlapConsequence->rank) >= 9));
		next if (defined $config->{novel_only} && $existing_vf ne '-');
		
		
		
		foreach my $tv(@{$new_vf->get_all_TranscriptVariations}) {
			
				# Conservation
				
				my %con = (
					  'type'		=> undef,
					  'gene_id'                  => undef,
					  'transcript_id'            => undef,
					  'protein_id'             => undef,
					  'uniprot_id'               => undef,
					  'uniprot_kb'			=> [],
					  'uniprot_trmbl'		=> [],
					  'omim'			=> undef,
					  'hgmd'			=> undef,
					  'entrez_gene_name'         => undef,
					  'entrez_gene_id'           => undef,
					  'entrez_gene_model_reactome' => undef,
					  'entrez_gene_model_summary' => undef,
					  'entrez_gene_model_other_source_anchor' => undef,
					  'entrez_gene_model_gene_rif' => undef,
					  'gene_description'         => undef,
					  'gene_ontology'            => undef,
					  'transcript_strand'        => $tv->transcript->strand(),
					  'transcript_SNP_Reference' => undef,
					  'transcript_SNP_Reads'     => undef,
					  'changes_protein'          => undef,
					  'amino_acid_reference'     => undef,
					  'amino_acid_variant'       => undef,
					  'altered_base_start'       => undef,
					  'altered_base_end'        => undef,
					  'altered_aa_start'         => undef,
					  'altered_aa_end'           => undef,
					  'protein_length'	=> undef,
					  'Model_Annotations'		=> [],
					  'Other_Annotations'		=>[],
					  'aligned_homolog_residues' => [],
					  'orthology'		=> [],
					  'reference_region'       => undef,
					  'protein_sequence'         => undef,
					  'homolog_species'          => []
					,    #species used for SNP residue conservation analysis
					  'alignment_score_change' => undef,    #see docs
					  'reference_amino_acid_conservation' =>
					    undef
					, #the percentage of ortholog residues identical to the reference
					 #equivalent to "Cindent" from Andreas Kowarsch et al, PLOS Computational Biology, 6:e1000923
				 	  'reference_amino_acid_score' =>
					    undef
					, #sum over scoring matrix scores for each pair between reference and aligned residues
					 #from homologs divided by the score you would get if all aligned residues match reference.
					 #equivalent to "Cblosum" from Andreas Kowarsch et al, PLOS Computational Biology, 6:e1000923
					  'context_conservation' =>
					    undef
					, #the average percent identity of protein region surrounding SNP, including SNP portion
					  'overlapping_protein_domains' =>
					    undef
					, #domains (from Ensembl) that overlap with SNP region on protein
					  'uniprot_overlapping_protein_features' => []
					, #protein features (from UniProt) that overlap with SNP region
					 'uniprot_parsed' => undef,
					  'uniprot_overlapping_protein_features_orthologues' => []
					, #protein features (from UniProt) on orthologous proteins that overlap with SNP region
					  'aligned_model_orthologue_proteins' =>undef
					,#holds hashes providing information on orthologous proteins from model
					  'disease' => undef,
					  #holds disease annotations information from uniprot_parsing
					model_phenotype_annotations => []
					, #whether changes in SNP region in model species are associated with a phenotype
					is_correlated_family =>undef,
					#whether the affected amino acid position is correlated with other amino acid
					#positions in protein family, according to the OMES procedure described in
					#Fodor and Aldrich PROTEINS: Structure, Function, and Bioinformatics 56:211Ð221\
					#Andreas Kowarsch et al, PLOS Computational Biology, 6:e1000923
					is_correlated_orthologues =>undef,
					#whether the affected amino acid position is correlated with other amino acid
					#positions in orthologues, according to the OMES procedure described in
					#Fodor and Aldrich PROTEINS: Structure, Function, and Bioinformatics 56:211Ð221 and
					#Andreas Kowarsch et al, PLOS Computational Biology, 6:e1000923
					model_orthologue_gene_id =>undef,    #the Ensembl Gene ID of the model orthologue
					entrez_gene_model_kegg_pathways =>undef,  #kegg pathways from Entrez Gene for model gene
					gerp_score => undef,
					entrez_gene_model_phenotypes =>undef    #phenotypes from Entrez Gene for model gene
					#entrez_gene_model_interactions_count =>undef #number of interactions from Entrez Gene for model gene
					#these are not just protein-protein interactions (could be protein binding promoter for example)
					);

			next if(defined $config->{coding_only} && !($tv->affects_transcript));
			#my @transcript_variation=$tv->get_all_alternate_TranscriptVariationAlleles;
			
			#foreach my $tva(@{$tv->get_all_alternate_TranscriptVariationAlleles}) {
				#my $type = join ",", map {$_->$method_name || $_->display_term} @{$tva->get_all_OverlapConsequences};
				my $method_name = $config->{terms}.'_term';
				my $type = join ",", @{$tv->consequence_type};
				$con{type}=$type;
				#my $is_damaging=undef;
				#if ( $type =~/SPLICE/ || $type =~/FRAMESHIFT/ || $type =='NON_SYNONYMOUS_CODING' || $type =~/STOP/ || $type == 'COMPLEX_INDEL' || $type =~ /UNKNOWN/) {
				#	$is_damaging=1;	
				#}
				#my $type;

				#my $gene = ($tv->transcript ? $config->{ga}->fetch_by_transcript_stable_id($tv->transcript->stable_id) : undef) unless defined $config->{whole_genome};
				my $gene = ($tv->transcript ? $config->{ga}->fetch_by_transcript_stable_id($tv->transcript->stable_id) : undef);
				my $extra;
				my $hgnc;
				
				
				
				
				# HGNC
				if(defined $config->{hgnc} && $gene) {
				  my @entries = grep {$_->database eq 'HGNC'} @{$gene->get_all_DBEntries()};
				  if(scalar @entries) {
					#$extra .= 'HGNC='.$entries[0]->display_id.';';
					$hgnc=$entries[0]->display_id;
				  }
				}
				
				# protein ID
				if(defined $config->{protein} && $tv->transcript->translation) {
					$extra .= 'ENSP='.$tv->transcript->translation->stable_id.';';
				}
				
				# HGVS
				#if(defined $config->{hgvs}) {
				#	$extra .= 'HGVSc='.$tva->hgvs_coding.';' if defined($tva->hgvs_coding);
				#	$extra .= 'HGVSp='.$tva->hgvs_protein.';' if defined($tva->hgvs_protein);
				#}
				
				# Conservation
				
				$con{transcript_strand} = $tv->transcript->strand();
				if ( defined( $tv->transcript->translation ) ) {
					$con{protein_id} = $tv->transcript->translation->stable_id();
				}
				$con{gene_id}          = $gene->stable_id();
				$con{gene_description} = $gene->description();
				#debug("before xrefs");
				
				my $xrefs = $tv->transcript->get_all_DBLinks();
				my @go_names = ();
				foreach my $xref ( @{$xrefs} ) {
					if ( $xref->dbname() eq 'EntrezGene' ) {
						$con{entrez_gene_name} = $xref->display_id();
						$con{entrez_gene_id}   = $xref->primary_id();
					}
					elsif ( $xref->dbname() eq 'Uniprot/SWISSPROT' ) {
				
						push(
						     @{ $con{uniprot_id} },
						     $xref->display_id()
						    );
					}
					elsif ( $xref->dbname() =~/^MIM/) {
						push(
						     @{ $con{omim} },
						     $xref->display_id()
						    );
					}
					
					elsif ( $xref->dbname() =~ /HGMD/) {
						push(
						     @{ $con{hgmd} },
						     $xref->display_id()
						    );
					}
						elsif ( $xref->dbname() =~ /ImmunoDB/) {
						push(
						     @{ $con{hgmd} },
						     $xref->display_id()
						    );
					}
					#elsif ( $xref->dbname() eq 'goslim_goa' ) {
					elsif ( $xref->dbname() =~ /GO/ ) {
						my $go_id = $xref->display_id();
						my $go_term;
						my $go_name;
						my $go_definition;
						
						if ( defined( $config->{goa} ) ) {
							$go_term = $config->{goa}->fetch_by_accession($go_id);
							$go_name       = $go_term->name();
							$go_definition = $go_term->definition();
						}
						if ( defined($go_name) ) {
							push( @go_names, "$go_name" );
						}
						else {
							push( @go_names, "$go_id" );
						}
					}
				}
				#debug("after xref");
				$con{gene_ontology}= join( ';', @{ get_unique( \@go_names ) } );
				if ( defined( $tv->pep_allele_string() ) ) {
					#reference is first
					if ( $tv->pep_allele_string =~ m/([^\/]+)\/([^\/]+)/ ) {
						$con{amino_acid_reference} = $1;
						$con{amino_acid_variant}   = $2;
						if ( $con{amino_acid_reference} ne $con{amino_acid_variant} ) {
							$con{changes_protein} = 1;
						}
					}
				}
				#debug("gene ontology");
				if ($con{changes_protein}) {
					if ( defined( $tv->transcript->translate() ) ) {
						$con{protein_sequence}= $tv->transcript->translate()->seq();
					}
					$con{altered_aa_start}   = $tv->translation_start();
					$con{altered_aa_end}     = $tv->translation_end();
					$con{altered_base_start} = $tv->cdna_start();
					$con{altered_base_end}   = $tv->cdna_end();
					$con{transcript_id} = $tv->transcript->stable_id;
					
					#$con{aligned_model_orthologue_proteins}= get_model_orthologous_residues($config, \%con);
					#debug("get model orthologue residues");
				}
				
				if ( $config->{use_uniprot} && $con{uniprot_id} ) {
					get_overlapping_protein_features_from_UniProt( $config, \%con, \%parsed );
					get_overlapping_protein_features_from_UniProt_orthologues($config, \%con, \%parsed);
				}
				#debug("uniprot");
				
				if ( $config->{use_ncbi} && $con{entrez_gene_id} ) {
					get_model_orthologue_info_from_NCBI( $config, \%con , \%parsed);
				}
				#debug("ncbi");
				
				
				my $scoring_matrix = File::Spec->catfile($ENV{NGS_SNP}, 'scripts', 'annotate_SNPs', 'data', 'blosum62.mat' );
				my $parser = Bio::Matrix::IO->new(
								  -format => 'scoring',
								  -file   => $scoring_matrix
								 );
				my $matrix = $parser->next_matrix;
				
				
				my $orthology;
				if ( ( $con{changes_protein} ) ) {
					$orthology = determine_aligned_homolog_residues( $config, \%con );
					#determine_alignment_score_change( $matrix, \%con );
					#determine_reference_amino_acid_conservation( \%con );
					#determine_reference_amino_acid_score( $matrix,\%con );
					##determine_context_average_percent_identity( $config,\%con );
				}
				
				if ( defined( $config->{model} ) ) {
					get_model_phenotypes( $config, \%con, $new_vf );
					#debug("get model  phenotypes ");
				}
				
				#$con{gerp_score}=get_conservation_score($config, \%con, $new_vf);
				
				my $model_an=determine_model_annotations( $config, \%con);
				my $other_an=determine_other_annotations( $config, \%con);
				
				#Information about protein site conservation
				
				my ($Amino_Acids_In_Orthologues, $Orthologue_Species, $Alignment_Score_Change, $C_blosum, $C_ident, $Model_Annotations, $GO, $interactions, $allele_frequencies, $allele_populations, $disease);
				
				#if (   ( defined( $con{aligned_homolog_residues} ) ) && ( scalar( @{ $con{aligned_homolog_residues} } ) > 0 ) ) {
				#	$Amino_Acids_In_Orthologues = join( '', @{ $con{aligned_homolog_residues} } );
				#}
				#
				#if (   ( defined( $con{homolog_species} ) ) && ( scalar( @{ $con{homolog_species} } ) > 0 ) ) {
				#	$Orthologue_Species = join( ';', @{ $con{homolog_species} } );
				#}
				
				
				
				if (   ( defined( $orthology ) ) ) {
					$Orthologue_Species = join( ';',  keys %$orthology );
					$Amino_Acids_In_Orthologues = join('',values %$orthology);	
				}
				
				if (   ( defined( %$all_freq ) ) && ( scalar( %$all_freq ) > 0 ) ) {
					$allele_populations = join(';', keys %$all_freq);
					$allele_frequencies = join(';', values %$all_freq);
					my @product;
					my $i = 0;
					if ($#freqs == $#counts) {
						  foreach(@freqs) {
						    $product[$i] = $_ * $counts[$i];
						    ++$i;
						  }
						}
				}
			
				###Add up allele frequencies
				
				#if (   ( defined( $con{entrez_gene_model_interactions_count} ) ) && ( scalar( @{ $con{entrez_gene_model_interactions_count} } ) > 0 ) ) {
					#$interactions = join( ';', @{ $con{entrez_gene_model_interactions_count} } );
				#}
				
				#$Alignment_Score_Change = $con{alignment_score_change};
				#$C_blosum = $con{reference_amino_acid_score};
				#$C_ident  = $con{reference_amino_acid_conservation};
				$GO = $con{gene_ontology};
				my $uniprot_annotations=$con{uniprot_parsed};
				if (   ( defined( $con{disease} ) ) && ( scalar( @{ $con{disease} } ) > 0 ) ) {
					$disease = join( ';', @{ $con{disease} } );
					}
				my @desc = split(/\[Source/, $con{gene_description});
				
				#  $Amino_Acids_In_Orthologues. "\n".
				#  $Orthologue_Species. "\n".
				#  $Alignment_Score_Change. "\n".
				#  $Context_Conservation. "\n".
				#  $C_blosum. "\n".
				#  $C_ident. "\n".
				#  $GO. "\n";
				#  #$Model_Annotations;
				foreach my $tva(@{$tv->get_all_alternate_TranscriptVariationAlleles}) {
					#$type = join ",", map {$_->$method_name || $_->display_term} @{$tva->get_all_OverlapConsequences};
					foreach my $tool (qw(SIFT PolyPhen Condel)) {
						my $lc_tool = lc($tool);
						
						if (my $opt = $config->{$lc_tool}) {
							my $want_pred   = $opt =~ /^p/i;
							my $want_score  = $opt =~ /^s/i;
							my $want_both   = $opt =~ /^b/i;
							
							if ($want_both) {
								$want_pred  = 1;
								$want_score = 1;
							}
							
							next unless $want_pred || $want_score;
							
							my $pred_meth   = $lc_tool.'_prediction';
							my $score_meth  = $lc_tool.'_score';
							
							my $pred = $tva->$pred_meth;
							
							if($pred) {
								$extra .= "$tool=";
								
								if ($want_pred) {
									$pred =~ s/\s+/\_/;
									$extra .= $pred;
								}
									
								if ($want_score) {
									my $score = $tva->$score_meth;
									
									if(defined $score) {
										if($want_pred) {
											$extra .= "($score)";
										}
										else {
											$extra .= $score;
										}
									}
								}
								
								$extra .= ';';
							}
						}
					}
					#debug("protein prediction tools ");
				
				
				$extra =~ s/\;$//g;
				
				#Genotypes & Ref/Alt_DP
				my $curidx=0;
				my @gts;
				my @ref_alts;
				my @abs;
			#	my $gt_column;
			#	my $refalt_column;
				my @info=split(",", $new_vf->source() );
				while ($curidx < $config->{number_of_samples}) {
					if ($info[$curidx] =~ /:/) {
						my ($each_gt, $each_dp, $ab) = split(":", $info[$curidx],3);
						push @gts, $each_gt;
						push @ref_alts, $each_dp;
						push @abs, $ab;
					}
					else {
						push @gts, $info[$curidx];
						push @ref_alts, ".";
						push @abs, ".";
					}
					$curidx++;
				}
				
				
				
				print $out_file_handle
				  $new_vf->variation_name, "\t",
				  $new_vf->seq_region_name, ":",
				  &format_coords($new_vf->start, $new_vf->end), "\t",
				  (split /\//, $new_vf->allele_string)[0], "\t",
				  (split /\//, $new_vf->allele_string)[1], "\t",
				  #$tva->variation_feature_seq, "\t",
				  #($gene ? $gene->stable_id : '.'), "\t",
				  #$new_vf->source, "\t",
				(scalar @gts ? join(",", @gts) : "."), "\t", 
				  (scalar @ref_alts ? join(",", @ref_alts) : "." ), "\t",
				  (scalar @abs ? join(",",@abs) : "."), "\t",
				  $filter, "\t",
				  $qual, "\t",
				  $type, "\t",
				  ($tv->transcript ? $tv->transcript->stable_id : "."), "\t",
				  ($hgnc || "."), "\t",
				  $desc[0], "\t",
				  $existing_vf, "\t",
				  #($allele_frequencies || '.'), "\t",
				  #($allele_populations || '.'), "\t",
				  ($alt_frequency || '.'), "\t",
				  &format_coords($tv->cdna_start, $tv->cdna_end), "\t",
				  &format_coords($tv->cds_start, $tv->cds_end), "\t",
				  &format_coords($tv->translation_start, $tv->translation_end), "\t",
				  ($tva->pep_allele_string || "."), "\t",
				  ($tva->display_codon_allele_string || "."), "\t",
				  ($extra || '.'), "\t",
				  ($Amino_Acids_In_Orthologues || '.'), "\t",
				  ($Orthologue_Species || '.'), "\t",
				  ($disease || '.'), "\t",
				  #($Alignment_Score_Change || '.'), "\t",
				  #($C_blosum || '.'), "\t",
				  #($C_ident || '.'), "\t",
				#  ($con{gerp_score} || '.'),"\t",
				  ($GO || '.'), "\t",
				  ($model_an || '.'), "\t",
				  ($other_an || '.'), "\t",
				  ($uniprot_annotations || '.'), "\n";
			
			}
		}
	}
	
}
sub configure {
	my $args = shift;
	
	my $config = {};
	
	GetOptions(
		$config,
		'help',
		
		# input options,
		'config=s',
		'input_file=s',
		'format=s',
		'number_of_samples=i',
		'use_uniprot',
		'use_ncbi',
		
		# DB options
		'species=s',
		'local',
		'model=s',
		'registry=s',
		'host=s',
		'user=s',
		'port=s',
		'password=s',
		'db_version=i',
		'genomes',
		
		# runtime options
		'most_severe',
		'buffer_size=i',
		'chunk_size=s',
		'check_ref',
		'check_existing=i',
		'failed=i',
		'whole_genome',
		'tmp_dir=s',
		'gp',
		
		# output options
		'output_file=s',
		'terms=s',
		'verbose',
		'quiet',
		'coding_only',
		'damaging',
		'novel_only',
		'protein',
		'hgnc',
		'hgvs',
		'sift=s',
		'polyphen=s',
		'condel=s',
		
		#Scripts
		'ncbi_script',
		'uniprot_script'
	);
	
	# print usage message if requested or no args supplied
	if(defined($config->{help}) || !$args) {
		&usage;
		exit(0);
	}
	
	# config file?
	if(defined $config->{config}) {
		
		open CONFIG, $config->{config} or die "ERROR: Could not open config file \"".$config->{config}."\"\n";
		
		while(<CONFIG>) {
			next if /^\#/;
			my ($key, $value) = split /\s+|\=/;
			$key =~ s/^\-//g;
			$config->{$key} = $value unless defined $config->{$key};
		}
		
		close CONFIG;
	}
	
	# check file format
	if(defined $config->{format}) {
		die "ERROR: Unrecognised input format specified \"".$config->{format}."\"\n" unless $config->{format} =~ /pileup|vcf|guess/i;
	}
	
	# output term
	if(defined $config->{terms}) {
		die "ERROR: Unrecognised consequence term type specified \"".$config->{terms}."\" - must be one of ensembl, so, ncbi\n" unless $config->{terms} =~ /ensembl|display|so|ncbi/i;
		if($config->{terms} =~ /ensembl|display/i) {
			$config->{terms} = 'display';
		}
		else {
			$config->{terms} = uc($config->{terms});
		}
	}
	
	# check nsSNP tools
	foreach my $tool(grep {defined $config->{lc($_)}} qw(SIFT PolyPhen Condel)) {
		die "ERROR: Unrecognised option for $tool \"", $config->{lc($tool)}, "\" - must be one of p (prediction), s (score) or b (both)\n" unless $config->{lc($tool)} =~ /^(s|p|b)/;
	}
	
	# summarise options if verbose
	if(defined $config->{verbose}) {
		my $header =<<INTRO;
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

version 2.0

By Will McLaren (wm2\@ebi.ac.uk)

Configuration options:

INTRO
		print $header;
		
		my $max_length = (sort {$a <=> $b} map {length($_)} keys %$config)[-1];
		
		foreach my $key(sort keys %$config) {
			print $key.(' ' x (($max_length - length($key)) + 4)).$config->{$key}."\n";
		}
		
		print "\n".("-" x 20)."\n\n";
	}
	
	# connection settings for Ensembl Genomes
	if($config->{genomes}) {
		$config->{host} ||= 'mysql.ebi.ac.uk';
		$config->{port} ||= 4157;
	}
	
	else {
		$config->{species} ||= "homo_sapiens";
		$config->{host}    ||= 'cortex.local';
		$config->{user}    ||= 'ensembl';
		$config->{password} ||='ensembl';
		
		
	}
	#
	#if($config->{local}) {
	#	$config->{host} ||= '127.0.0.1';
	#	$config->{port} ||= '3308';
	#	$config->{user} ||= 'tpaige';
	#	$config->{password} ||= '';
	#}
	#else {
	#	$config->{host} ||= 'cortex.local';
	#	$config->{port} ||= '3306';
	#	$config->{user} ||= 'tpaige';
	#	$config->{password} ||= '';
	#}
	
	# connection settings for main Ensembl

	
	# set defaults
	$config->{user}         ||= 'ensembl';
	$config->{password}         ||= 'ensembl';
	$config->{buffer_size}  ||= 5000;
	$config->{chunk_size}   ||= '50kb';
	$config->{output_file}  ||= "variant_effect_output.txt";
	$config->{tmpdir}       ||= '/tmp';
	$config->{format}       ||= 'guess';
	$config->{terms}        ||= 'display';
	
	$config->{sift} = 'b' unless defined $config->{sift};
	$config->{polyphen} = 'b' unless defined $config->{polyphen};
	$config->{condel} = 'b' unless defined $config->{condel};
	$config->{hgnc} = 1;
	$config->{model} ||= 'homo_sapiens';
	$config->{include_failed} = 0 unless defined $config->{include_failed};
	$config->{check_existing} = 1 unless (defined $config->{check_existing} || defined $config->{whole_genome});
	$config->{chunk_size} =~ s/mb?/000000/i;
	$config->{chunk_size} =~ s/kb?/000/i;
	
	#scripts
	$config->{ncbi_script} = File::Spec->catfile(
        $ENV{NGS_SNP},   'scripts',
        'annotate_SNPs', 'scripts',
        'ncbi_search',   'ncbi_search.pl'
	);
	
	$config->{uniprot_script} = File::Spec->catfile(
        $ENV{NGS_SNP}, 'scripts', 'annotate_SNPs', 'scripts',
        'ebi_fetch',   'ebi_fetch.pl'
	);
	# connect to databases
	$config->{reg} = &connect_to_dbs($config);
	
	
	# get input file handle
	$config->{in_file_handle} = &get_in_file_handle($config);
	
	# configure output file
	$config->{out_file_handle} = &get_out_file_handle($config);
	
	return $config;
}

sub usage {
	my $usage =<<END;
#----------------------------------#
# ENSEMBL VARIANT EFFECT PREDICTOR #
#----------------------------------#

version 2.0

By Will McLaren (wm2\@ebi.ac.uk)

Usage:
perl variant_effect_predictor.pl [arguments]

Options
--help                 Display this message and quit
--verbose              Display verbose output as the script runs [default: off]
--quiet                Suppress status and warning messages [default: off]

--config               Load configuration from file. Any command line options
                       specified overwrite those in the file [default: off]

-i | --input_file      Input file - if not specified, reads from STDIN. Files
                       may be gzip compressed.
--format               Alternative input file format - one of "pileup", "vcf"
-o | --output_file     Output file [default: "variant_effect_output.txt"]

-t | --terms           Type of consequence terms to output - one of "ensembl", "SO",
                       "NCBI" [default: ensembl]
					   
--sift=[p|s|b]         Add SIFT [p]rediction, [s]core or [b]oth [default: off]
--polyphen=[p|s|b]     Add PolyPhen [p]rediction, [s]core or [b]oth [default: on]
--condel=[p|s|b]       Add Condel SIFT/PolyPhen consensus [p]rediction, [s]core or
                       [b]oth [default: on]

NB: SIFT, PolyPhen and Condel predictions are currently available for human only

--number_of_samples	VCF file only.  Specifies the number of samples for which genotypes are present in an input VCF file.  Default is to guess.
--hgvs                 Output HGVS identifiers (coding and protein) [default: off]
--protein              Output Ensembl protein identifer [default: off]

--coding_only          Only return consequences that fall in the coding region of
                       transcripts (does not include SPLICE variants) [default: off]

--damaging 	Only return consequences that are most severe (splice sites, stop lost/gained, complex indels, non-synonymous coding, frameshift coding). [default: off]
		       
--novel_only		Only return consequences that do not have an existing rs number.  Disabled if run in wgs mode [default: off]	

--check_ref            If specified, checks supplied reference allele against stored
                       entry in Ensembl Core database [default: off]
--check_existing=[0|1] If set to 1, checks for existing co-located variations in the
                       Ensembl Variation database [default: 1]
--failed=[0|1]         If set to 1, includes variations flagged as failed when checking
                       for co-located variations. Only applies in Ensembl 61 or later
                       [default: 0]
--gp                   If specified, tries to read GRCh37 position from GP field in the
                       INFO column of a VCF file. Only applies when VCF is the input
                       format and human is the species [default: off]

-s | --species         Species to use [default: "human"]
-model [STRING]        model species to use when filling in the 
                      'Model_Annotations' column in the output. The model
                      species can be the same as the source species (Optional;
                      default value is 'homo_sapiens').

--host                 Manually define database host [default: "ensembldb.ensembl.org"]
-u | --user            Database username [default: "anonymous"]
--port                 Database port [default: 5306]
--password             Database password [default: no password]
--genomes              Sets DB connection params for Ensembl Genomes [default: off]
-r | --registry_file   Registry file to use defines DB connections [default: off]
                       Defining a registry file overrides above connection settings.
--db_version=[number]  Force script to load DBs from a specific Ensembl version. Not
                       advised due to likely incompatibilities between API and DB
					   
-w | --whole_genome    Run in whole genome mode [default: off]
                       Recommended for use with data covering a whole genome,
                       chromosome or gene/set of genes e.g. from resequencing. For
                       better performance, set --buffer_size higher (>10000 if memory
                       allows).
-b | --buffer_size     Sets the number of variants sent in each batch [default: 5000]
                       Increasing buffer size can retrieve results more quickly
                       but requires more memory. Only applies to whole genome mode.
--chunk_size           Sets the chunk size of internal data structure [default: 50kb]
                       Setting this lower may improve speed for variant-dense
                       datasets. Only applies to whole genome mode.
					   
NB: Whole genome mode disables --check_existing option and gene column by default.
END

	print $usage;
}


sub connect_to_dbs {
	my $config = shift;
	
	# get registry
	my $reg = 'Bio::EnsEMBL::Registry';
	
	# load DB options from registry file if given
	if(defined($config->{registry})) {
		debug("Loading DB config from registry file ", $config->{registry});
		$reg->load_all($config->{registry});
	}
	
	# otherwise manually connect to DB server
	else {
		#$reg->load_registry_from_db(
		#	-host       => '127.0.0.1',
		#	-user       => 'tpaige',
		#	-pass       => '',
		#	-port       => 3308,
		#	-db_version => 61,
		#	#-species    => $config->{species} =~ /^[a-z]+\_[a-z]+/i ? $config->{species} : undef,
		#	-verbose    => 1,
		#);	
		$reg->load_registry_from_db(
			-host       => $config->{host},
			-user       => $config->{user},
			-pass       => $config->{password},
			-port       => $config->{port},
			-db_version => $config->{db_version},
			-species    => $config->{species} =~ /^[a-z]+\_[a-z]+/i ? $config->{species} : undef,
			-verbose    => 1,
			-wait_timeout => 28800
		);
	}
	
	#$reg->set_disconnect_when_inactive();
	
	if($config->{verbose}) {
		# get a meta container adaptors to check version
		my $core_mca = $reg->get_adaptor($config->{species}, 'core', 'metacontainer');
		my $var_mca = $reg->get_adaptor($config->{species}, 'variation', 'metacontainer');
		
		if($core_mca && $var_mca) {
			debug(
				"Connected to core version ", $core_mca->get_schema_version, " database ",
				"and variation version ", $var_mca->get_schema_version, " database"
			);
		}
	}
	
	return $reg;
}


sub get_in_file_handle {
	my $config = shift;

	# define the filehandle to read input from
	my $in_file_handle = new FileHandle;
	
	if(defined($config->{input_file})) {
		
		# check defined input file exists
		die("ERROR: Could not find input file ", $config->{input_file}, "\n") unless -e $config->{input_file};
		
		if($config->{input_file} =~ /\.gz$/){
			$in_file_handle->open("zcat ". $config->{input_file} . " | " ) or die("ERROR: Could not read from input file ", $config->{input_file}, "\n");
		}
		else {
			$in_file_handle->open( $config->{input_file} ) or die("ERROR: Could not read from input file ", $config->{in_file}, "\n");
		}
	}
	
	# no file specified - try to read data off command line
	else {
		$in_file_handle = 'STDIN';
		debug("Reading input from STDIN (or maybe you forgot to specify an input file?)...") unless defined $config->{quiet};
	}
	
	return $in_file_handle;
}


sub get_out_file_handle {
	my $config = shift;
	
	# define filehandle to write to
	my $out_file_handle = new FileHandle;
	$out_file_handle->open(">".$config->{output_file}) or die("ERROR: Could not write to output file ", $config->{output_file}, "\n");
	
	# make header
	my $time = &get_time;
	my $core_mca = $config->{reg}->get_adaptor($config->{species}, 'core', 'metacontainer');
	my $db_string = $core_mca->dbc->dbname." on ".$core_mca->dbc->host if defined $core_mca;
	my $version_string =
		"Using API version ".$config->{reg}->software_version.
		", DB version ".(defined $core_mca && $core_mca->get_schema_version ? $core_mca->get_schema_version : '?');
	
	my $header =<<HEAD;
## ENSEMBL VARIANT EFFECT PREDICTOR v2.0
## Output produced at $time
## Connected to $db_string
## $version_string
## Extra column keys:
## HGNC     : HGNC gene identifier
## ENSP     : Ensembl protein identifer
## HGVSc    : HGVS coding sequence name
## HGVSp    : HGVS protein sequence name
## SIFT     : SIFT prediction
## PolyPhen : PolyPhen prediction
## Condel   : Condel SIFT/PolyPhen consensus prediction
HEAD
	
	# add headers
	print $out_file_handle $header;
	
	# add column headers
	print $out_file_handle join "\t", qw(
		#Uploaded_variation
		Location
		Ref
		Alt	
		Genotypes
		Ref/Alt_DP
		AB
		Filter
		Qual
		Consequence
		Transcript
		HGNC
		Description
		Existing_variation
		Allele_pop
		Average_Alt_Allele_Freq
		cDNA_position
		CDS_position
		Protein_position
		Amino_acids
		Codons
		Prediction
		AAs_in_Orthologues
		Orthologues
		Disease
		Ontology
		Model_Annotations
		Other_Annotations
		Uniprot_Annotations
	);
	
	print $out_file_handle "\n";
	
	return $out_file_handle;
}


sub parse_line {
	my $config = shift;
	my $line   = shift;
	
	my @data = (split /\s+/, $_);
	
	
	
	# pileup: chr1 60 T A
	if(
	   ($config->{format} =~ /pileup/i) ||
	   (
			$data[0] =~ /(chr)?\w+/ &&
			$data[1] =~ /\d+/ &&
			$data[2] =~ /^[ACGTN-]+$/ &&
			$data[3] =~ /^[ACGTNRYSWKM*+\/-]+$/
		)
	) {
		my @return = ();
		
		if($data[2] ne "*"){
			my $var;
			
			if($data[2] =~ /^[A|C|G|T]$/) {
				$var = $data[2];
			}
			else {
				($var = unambiguity_code($data[3])) =~ s/$data[2]//ig;
			}
			if(length($var)==1){
				push @return, [$data[0], $data[1], $data[1], $data[2]."/".$var, 1, undef];
			}
			else{
				for my $nt(split //,$var){
					push @return, [$data[0], $data[1], $data[1], $data[2]."/".$nt, 1, undef];
				}
			}
		}
		else{ #indel
			my @genotype=split /\//,$data[3];
			foreach my $allele(@genotype){
				if(substr($allele,0,1) eq "+") { #ins
					push @return, [$data[0], $data[1]+1, $data[1], "-/".substr($allele,1), 1, undef];
				}
				elsif(substr($allele,0,1) eq "-"){ #del
					push @return, [$data[0], $data[1], $data[1]+length($data[3])-4, substr($allele,1)."/-", 1, undef];			
				}
				elsif($allele ne "*"){
					warn("WARNING: invalid pileup indel genotype: $line\n") unless defined $config->{quiet};
					push @return, ['non-variant'];
				}
			}
		}
		return \@return;
	}
	
	# VCF: 20      14370   rs6054257 G     A      29    0       NS=58;DP=258;AF=0.786;DB;H2          GT:GQ:DP:HQ
	elsif(
		($config->{format} =~ /vcf/i) ||
		(
			$data[0] =~ /(chr)?\w+/ &&
			$data[1] =~ /\d+/ &&
			$data[3] =~ /^[ACGTN-]+$/ &&
			$data[4] =~ /^([\.ACGTN-]+\,?)+$/
		)
	) {
		$config->{format}='vcf';
		# non-variant line in VCF, return dummy line
		if($data[4] eq '.') {
			return [['non-variant']];
		}
		
		# get relevant data
		my ($chr, $start, $end, $ref, $alt) = ($data[0], $data[1], $data[1], $data[3], $data[4]);
		
		#Paige 5-1-11
		my $number_of_columns=scalar( @data );
		if ( !defined $config->{number_of_samples}) {
			$config->{number_of_samples}=$number_of_columns - 9;
		}
		
		## Get Genotypes
		my @format=(split ':', $data[8]);
		my $curidx = 0;
		my @genotypes;
		while ($curidx < $config->{number_of_samples}) {
			$genotypes[$curidx]=$data[(9+$curidx)];
			#my @info= (split ':', $genotypes[$curidx]);
			#my $g =$info[0];
			#my (@depth, $ref_dp, $alt_dp, $ab, $filt_depth, $gq, $pl, $tmp);
			my ($gt,$ab,$ad,$dp,$gq,$pl);
			if ($genotypes[$curidx] eq './.' ) {
				$gt="./."
			}
			elsif ($format[1] eq 'AD') {
				#unless ($info[0] eq "./.") {
				#unless ($genotypes[$curidx] == "./.") {
					#($ref_dp, $alt_dp) = (split ',', $info[1]);
					($gt,$ad,$dp,$gq,$pl)=split(":", $genotypes[$curidx]);
					$ab=".";
				}
			
			elsif ($format[2] eq 'AD') {
				#unless ($info[0] eq "./."){
				#unless ($genotypes[$curidx] == "./.") {
					#($ref_dp, $alt_dp) = (split ',', $info[2]);
					($gt,$ab,$ad,$dp,$gq,$pl)=split(":", $genotypes[$curidx]);
				}
			
			else {
				my @other=split /\:/, $genotypes[$curidx];
				$gt=$other[0];
			}
			#if (defined $ref_dp && defined $alt_dp) {
			if (defined $ad) {
				#$genotypes[$curidx]=$info[0].':'.$ref_dp.':'.$alt_dp;}
				$genotypes[$curidx]=$gt.":".join("/",split(",", $ad) ).":".$ab }
			else {
				$genotypes[$curidx]=$gt;
			}
			
			$curidx++;
		}
		 
		 
		my $gt = join(",", @genotypes);
		
		if(defined $config->{use_gp}) {
			$chr = undef;
			$start = undef;
			
			foreach my $pair(split /\;/, $data[7]) {
				my ($key, $value) = split /\=/, $pair;
				if($key eq 'GP') {
					($chr,$start) = split /\:/, $value;
					$end = $start;
				}
			}
			
			unless(defined($chr) and defined($start)) {
				warn "No GP flag found in INFO column" unless defined $config->{quiet};
				return [['non-variant']];
			}
		}
		
		# adjust end coord
		$end += (length($ref) - 1);
		
		# find out if any of the alt alleles make this an insertion or a deletion
		my ($is_indel, $is_sub, $ins_count, $total_count);
		foreach my $alt_allele(split /\,/, $alt) {
			$config->{is_indel} = 1 if $alt_allele =~ /D|I/;
			$is_indel = 1 if length($alt_allele) != length($ref);
			$is_sub = 1 if length($alt_allele) == length($ref);
			$ins_count++ if length($alt_allele) > length($ref);
			$total_count++;
		}
		
		# multiple alt alleles?
		if($alt =~ /\,/) {
			if($is_indel) {
				
				my @alts;
				
				if($alt =~ /D|I/) {
					foreach my $alt_allele(split /\,/, $alt) {
						# deletion (VCF <4)
						if($alt_allele =~ /D/) {
							push @alts, '-';
						}
						
						elsif($alt_allele =~ /I/) {
							$alt_allele =~ s/^I//g;
							push @alts, $alt_allele;
						}
					}
				}
				
				else {
					$ref = substr($ref, 1);
					$ref = '-' if $ref eq '';
					$start++;
					
					foreach my $alt_allele(split /\,/, $alt) {
						$alt_allele = substr($alt_allele, 1);
						$alt_allele = '-' if $alt_allele eq '';
						push @alts, $alt_allele;
					}
				}
				
				$alt = join "/", @alts;
			}
			
			else {
				# for substitutions we just need to replace ',' with '/' in $alt
				$alt =~ s/\,/\//;
			}
		}
		
		else {
			if($is_indel) {
				# deletion (VCF <4)
				if($alt =~ /D/) {
					my $num_deleted = $alt;
					$num_deleted =~ s/\D+//g;
					$end += $num_deleted - 1;
					$alt = "-";
					$ref .= ("N" x ($num_deleted - 1)) unless length($ref) > 1;
				}
				
				# insertion (VCF <4)
				elsif($alt =~ /I/) {
					$ref = '-';
					$alt =~ s/^I//g;
					$start++;
				}
				
				# insertion or deletion (VCF 4+)
				else {
					# chop off first base
					$ref = substr($ref, 1);
					$alt = substr($alt, 1);
					
					$start++;
					
					if($ref eq '') {
						# make ref '-' if no ref allele left
						$ref = '-';
					}
					
					# make alt '-' if no alt allele left
					$alt = '-' if $alt eq '';
				}
			}
		}
		
		return [[$chr, $start, $end, $ref."/".$alt, 1, ($data[2] eq '.' ? undef : $data[2]), $gt]];
		
	}
	
	# our format
	else {
		# we allow commas as delimiter so re-split
		@data = (split /\s+|\,/, $_);
		return [\@data];
	}
}


sub whole_genome_fetch {
	my $config = shift;
	my $vf_hash = shift;
	my $transcript_cache = shift;
	
	my $up_down_size = MAX_DISTANCE_FROM_TRANSCRIPT;
	
	my %vf_done;
	my @return;
	
	foreach my $chr(keys %$vf_hash) {
		my $slice;
		
 		# first try to get a chromosome
 		eval { $slice = $config->{sa}->fetch_by_region('chromosome', $chr); };
 		
 		# if failed, try to get any seq region
 		if(!defined($slice)) {
 			$slice = $config->{sa}->fetch_by_region(undef, $chr);
 		}
		
		debug("Analyzing chromosome $chr") unless defined($config->{quiet});
		
		$transcript_cache->{$chr} = $slice->get_all_Transcripts() unless defined($transcript_cache->{$chr});
		
		my $tr_count = scalar @{$transcript_cache->{$chr}};
		
		debug("Fetched $tr_count transcripts") unless defined($config->{quiet});
		
		my $tr_counter;
		
		while($tr_counter < $tr_count) {
			
			my $tr = $transcript_cache->{$chr}->[$tr_counter++];
			
			if($tr_counter =~ /000$/) {
				#debug("Analysed $tr_counter\/$tr_count transcripts") unless defined($config->{quiet});
			}
			
			# do each overlapping VF
			my $chr = $tr->seq_region_name;
			my $s = $tr->start - $up_down_size;
			my $e = $tr->end + $up_down_size;
			
			# get the chunks this transcript overlaps
			my %chunks;
			$chunks{$_} = 1 for (int($s/$config->{chunk_size})..int($e/$config->{chunk_size}));
			map {delete $chunks{$_} unless defined($vf_hash->{$chr}{$_})} keys %chunks;
			
			foreach my $chunk(keys %chunks) {				
				foreach my $pos(grep {$_ >= $s && $_ <= $e} keys %{$vf_hash->{$chr}{$chunk}}) {
					foreach my $vf(@{$vf_hash->{$chr}{$chunk}{$pos}}) {
						$vf->add_TranscriptVariation(Bio::EnsEMBL::Variation::TranscriptVariation->new(
							-transcript => $tr,
							-variation_feature => $vf,
							-adaptor => $config->{tva},
						));
						
						if(!defined($vf_done{$vf->{'variation_name'}.'_'.$vf->{'start'}})) {
							push @return, $vf;
							$vf_done{$vf->{'variation_name'}.'_'.$vf->{'start'}} = 1;
						}
					}
				}
			}
		}
		
		# clean hash
		delete $vf_hash->{$chr};
	}
	
	debug("Calculating and writing output") unless defined($config->{quiet});
	
	
	return \@return;
}

# gets time
sub get_time() {
	my @time = localtime(time());

	# increment the month (Jan = 0)
	$time[4]++;

	# add leading zeroes as required
	for my $i(0..4) {
		$time[$i] = "0".$time[$i] if $time[$i] < 10;
	}

	# put the components together in a string
	my $time =
 		($time[5] + 1900)."-".
 		$time[4]."-".
 		$time[3]." ".
		$time[2].":".
		$time[1].":".
		$time[0];

	return $time;
}



# prints debug output with time
sub debug {

	my $text = (@_ ? (join "", @_) : "No message");
	my $time = get_time;

	print $time." - ".$text.($text =~ /\n$/ ? "" : "\n");
}

sub format_coords {
	my ($start, $end) = @_;
	
	if(!defined($start)) {
		return '.';
	}
	elsif(!defined($end)) {
		return $start;
	}
	elsif($start == $end) {
		return $start;
	}
	elsif($start > $end) {
		return $end.'-'.$start;
	}
	else {
		return $start.'-'.$end;
	}
}

sub find_existing {
	my $new_vf = shift;

	if(defined($new_vf->adaptor->db)) {
		
		my $sth1 = $new_vf->adaptor->db->dbc->prepare(qq{
		  SELECT variation_name, source_id
		  FROM variation_feature
		  WHERE seq_region_id = ?
		  AND seq_region_start = ?
		  AND seq_region_end = ?
		});
		
	
		

		$sth1->execute($new_vf->slice->get_seq_region_id, $new_vf->start, $new_vf->end);
		#$sth2->execute($new_vf->slice->get_seq_region_id, $new_vf->start, $new_vf->end);
		
		#my $vid=$new_vf->{dbID};
		#print $new_vf->display_id;
		#print $new_vf->id;
		#my $var_name=$new_vf->variation_name;
		#my $new_var=$new_vf->variation();
		#
		#my $var = $config->{va}->fetch_by_name('rs79585140');
		#
		my ($varia_name, $source, $allele, $count, $pop_id, $sample_name, $rs, $name);
		$sth1->bind_columns( \$varia_name, \$source );
		#$sth2->bind_columns( \$allele, \$count, \$pop_id, \$sample_name, \$rs );

		#
		#my ($name, $source);
		#$sth->bind_columns(\$name, \$source);
		
		my %by_source;
		my %by_source2;

		my @frequencies;
		
		push @{$by_source{$source}}, $varia_name while $sth1->fetch;
		$sth1->finish;
		
		#push @{$by_source2{$allele}}, $count while $sth2->fetch;
		#$sth2->finish;
		
		if(scalar keys %by_source) {
			foreach my $s(sort {$a <=> $b} keys %by_source) {
				return shift @{$by_source{$s}};
			}
		}
	}
	
	return undef;
}

#sub get_model_orthologous_residues {
#    my $config = shift;
#    my $con = shift;
#
#    if (   ( !defined( $con->{gene_id} ) )
#        || ( !defined( $con->{protein_id} ) )
#        || ( !defined( $con->{altered_aa_start} ) )
#        || ( !defined( $config->{translation_adaptor} ) ) )
#    {
#        return;
#    }
#    
##    my $comparadb= new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
##	-host => 'ensembldb.ensembl.org',
##	-port => 5306,
##	-user => 'anonymous',
##	-dbname => 'ensembl_compara_62',
##    );
#    my @model_orthologous_residues = ();
#
#    #$options->{ma} is a Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor
#    #$member is a Bio::EnsEMBL::Compara::Member
#
#    my $member = $config->{ma}->fetch_by_source_stable_id( 'ENSEMBLGENE',$con->{gene_id}  );
#
#    #Rarely the $member object may be undef
#    if ( !defined($member) ) {
#        return;
#    }
#
#    #$options->{ha} is a Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor
#    #$homologies is a list of Bio::EnsEMBL::Compara::Homology objects
#
#	my $homologies = [];
# 
#	#my $homologies_tut=$config->{ha}->fetch_all_by_Member($member);
#	
#
#    ##my $all_genome_dbs = $config->{genome_db_adaptor}->fetch_all();
#    ##my $fetch_all_mlss = $config->{mlss_adaptor}->fetch_all;
#    ###my $gdb = $config->{genome_db_adaptor}->fetch_by_name_assembly("Homo sapiens", 'NCBI36');
#    ###my $method_link_species_sets = $config->{mlss_adaptor}->fetch_all_by_genome_db($gdb);
#    ##my $method_link_species_set =$config->{mlss_adaptor}->fetch_by_method_link_type_species_set_name("EPO", "mammals");
#    ##my $gerp_const_mlss = $config->{mlss_adaptor}->fetch_by_method_link_type_species_set_name("GERP_CONSTRAINED_ELEMENT", "mammals");
#    ##my $gerp_conserv_mlss = $config->{mlss_adaptor}->fetch_by_method_link_type_species_set_name("GERP_CONSERVATION_SCORE", "mammals");
#   
#
#    push(
#        @{$homologies},
#        @{  $config->{ha}->fetch_all_by_Member_paired_species( $member,
#                'homo_sapiens', ['ENSEMBL_ORTHOLOGUES'] )
#            }
#    );
#    
#    foreach my $homology ( @{$homologies} ) {
#	
#	
#
#     #$homologues is an array ref of (2) Bio::EnsEMBL::Compara::Member objects
#     #$aln is Bio::SimpleAlign object
#	
#
#        #added eval block on 2009-10-22 because
#        #$aln = $homology->get_SimpleAlign();
#        #throws exception for ENSBTAT00000060548
#        my $homologues = undef;
#        my $taxon1     = undef;
#        my $taxon2     = undef;
#        my $aln        = undef;
#
#        eval {
#            $homologues = $homology->gene_list();
#            $taxon1     = $$homologues[0]->taxon;
#            $taxon2     = $$homologues[1]->taxon;
#            $aln        = $homology->get_SimpleAlign();
#        };
#        if ($@) {
#            next;
#        }
#
#        if ( !( defined($aln) ) ) {
#            next;
#        }
#
#        #confirm that reference protein is in the alignment
#        if (!(  scalar( $aln->each_seq_with_id( $con->{protein_id} )) == 1 ))
#        {
#            next;
#        }
#
#    #get the column containing the residue encoded by the SNP-containing codon
#        my $col = undef;
#
#     #column_from_residue_number can throw an exception if the altered residue
#     #in the reference is the one that encodes the stop codon
#        eval {
#            $col = $aln->column_from_residue_number(
#                $con->{protein_id},
#                $con->{altered_aa_start}
#            );
#        };
#        if ($@) {
#            next;
#        }
#
#        #get an alignment containing only the desired column
#        #my $sub_align = $aln->slice( $col, $col );
#        my $sub_align = undef;
#        eval { $sub_align = $aln->slice( $col, $col ); };
#        if ($@) {
#            next;
#        }
#
##confirm that there are two sequences in the subalignment
##(sequences are excluded if the alignment slice  from the sequence)
#        if ( !( scalar( $sub_align->each_seq() ) == 2 ) ) {
#            next;
#        }
#
#        my $seq1 = $sub_align->get_seq_by_pos(1);
#        my $seq2 = $sub_align->get_seq_by_pos(2);
#
#        #in all the examples I've examined, $seq1 is the reference sequence
#        if ( $seq1->display_id eq $con->{protein_id} ) {
#		
#	
#		my $model_gene = $$homologues[1]->get_Gene;
#		my $entrez_gene_id;
#		my $refseq_protein_id;
#		foreach my $xref ( @{ $model_gene->get_all_DBLinks() } ) {
#			if ( $xref->dbname() eq 'EntrezGene' ) {
#				$entrez_gene_id = $xref->primary_id();
#			}
#			elsif ( $xref->dbname() eq 'RefSeq_peptide' ) {
#				$refseq_protein_id = $xref->display_id();
#			}
#		}
#		#obtain the the Ensembl ID of the model protein
#		my $model_translation_ID = $seq2->display_id;
#		#obtain the position in the model protein of the amino acid that aligns with the SNP-altered residue
#		my $model_aligned_aa_start = $seq2->start;
#		my $model_aligned_aa_end   = $seq2->end;
#		
#		my %orthologue = (
#                id             => $model_translation_ID,
#                start          => $model_aligned_aa_start,
#                end            => $model_aligned_aa_end,
#                sequence       => undef,
#                entrez_gene_id => $entrez_gene_id,
#                refseq_peptide =>
#                    undef,   #will only be used to try to get NCBI Gene record
#                uniprot_id => [],
#                go         => undef
#		);
#		
#		#attempt to fill in usefull IDs and other info
#		my $model_translation_object = $config->{translation_adaptor}->fetch_by_stable_id($model_translation_ID);
#		$orthologue{sequence} = $model_translation_object->seq();
#		
#		my $xrefs = $model_translation_object->get_all_DBLinks();
#		
#		my @go_names = ();
#		foreach my $xref ( @{$xrefs} ) {
#			if ( $xref->dbname() eq 'Uniprot/SWISSPROT' ) {
#				push( @{ $orthologue{uniprot_id} }, $xref->display_id() );
#				}
#			elsif ( $xref->dbname() eq 'RefSeq_peptide' ) {
#				$orthologue{refseq_peptide} = $xref->display_id();
#			}
#			elsif ( $xref->dbname() eq 'goslim_goa' ) {
#				my $go_id = $xref->display_id();
#				my $go_term;
#				my $go_name;
#				my $go_definition;
#				#if ( defined( $options->{goa} ) ) {
#				$go_term = $config->{goa}->fetch_by_accession($go_id);
#				$go_name       = $go_term->name();
#				$go_definition = $go_term->definition();
#				#}
#				if ( defined($go_name) ) {
#					push( @go_names, "[$go_id]:$go_name" );
#				}
#				else {
#					push( @go_names, "$go_id" );
#				}
#			}
#		}
#		$orthologue{go} = get_unique( \@go_names );
#		if ( !defined( $orthologue{refseq_peptide} ) ) {
#			$orthologue{refseq_peptide} = $refseq_protein_id;
#		}
#		
#		push( @model_orthologous_residues, \%orthologue );
#	}
#	else {
#            die("Unexpected sequence ordering in alignment.");
#	}
#    }
#    return \@model_orthologous_residues;
#}
	

sub determine_aligned_homolog_residues {
    my $config = shift;
    my $con = shift;
    
    if ( !( $con->{changes_protein} ) ) {
        return;
    }
    
    #$options->{ma} is a Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor
    #$member is a Bio::EnsEMBL::Compara::Member
    my $member = $config->{ma}->fetch_by_source_stable_id( 'ENSEMBLGENE', $con->{gene_id} );

    #Rarely the $member object may be undef
    if ( !defined($member) ) {
        return;
    }

    #$options->{ha} is a Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor
    #$homologies is a list of Bio::EnsEMBL::Compara::Homology objects
    my $homologies = [];
    
    #obtain orthologues from all species
    push(
	 @{$homologies},
	 @{   $config->{ha}->fetch_all_by_Member_method_link_type( $member, 'ENSEMBL_ORTHOLOGUES' )
	  }
	);
	
	my %results;
#	    my %results= (
#		species                            => [],
#		residue                      => []
#    );
#	    

    
    
    foreach my $homology ( @{$homologies} ) {
	
	#$homology->print_homology;
	

     #$homologues is an array ref of (2) Bio::EnsEMBL::Compara::Member objects
     #$aln is Bio::SimpleAlign object

        #added eval block on 2009-10-22 because
        #$aln = $homology->get_SimpleAlign();
        #throws exception for ENSBTAT00000060548
        my $homologues = undef;
        my $taxon1     = undef;
        my $taxon2     = undef;
        my $aln        = undef;



        eval {
            $homologues = $homology->gene_list();
            $taxon1     = $$homologues[0]->taxon;
            $taxon2     = $$homologues[1]->taxon;
            $aln        = $homology->get_SimpleAlign();
        };
        if ($@) {
            next;
        }

        if ( !( defined($aln) ) ) {
            next;
        }

        #confirm that reference protein is in the alignment

        if (!(  scalar(
                    $aln->each_seq_with_id($con->{protein_id} )
                ) == 1
            )
            )
        {
            next;
        }

    #get the column containing the residue encoded by the SNP-containing codon
        my $col_start = undef;
        my $col_end = undef;

    eval {
            $col_start = $aln->column_from_residue_number(
                $con->{protein_id},
                $con->{altered_aa_start}
            );
            $col_end = $aln->column_from_residue_number(
                $con->{protein_id},
                $con->{altered_aa_end}
                )
        };
        if ($@) {
            next;
        }

        #get an alignment containing only the desired column
        #my $sub_align = $aln->slice( $col, $col );
        my $sub_align = undef;
        eval { $sub_align = $aln->slice( $col_start, $col_end ); };
        if ($@) {
            next;
        }
        #confirm that there are two sequences in the subalignment
#(sequences are excluded if the alignment slice contains no residues from the sequence)
        if ( !( scalar( $sub_align->each_seq() ) == 2 ) ) {
            next;
        }

        my $seq1 = $sub_align->get_seq_by_pos(1);
        my $seq2 = $sub_align->get_seq_by_pos(2);

        #in all the examples I've examined, $seq1 is the reference sequence
        if ( $seq1->display_id eq $con->{protein_id} ) {
            if (uc( $seq1->seq() ) ne
                uc( $con->{amino_acid_reference} ) )
            {
                my $transpose=uc( $seq1->seq() );
                $transpose =~ tr/-//d;
                if ($transpose ne uc( $con->{amino_acid_reference} ) ) {
                    
    #if the residue in the alignment is selenocysteine (U), change the reference to U
    #as U in reference seems to be given as * in string returned by $con->pep_allele_string()
                    if (   ( uc( $seq1->seq() ) eq 'U' )
                        && ( $con->{amino_acid_reference} eq '*' ) )
                    {
                        $con->{amino_acid_reference} = 'U';
    
              #assume base change in U-coding residue is NON_SYNONYMOUS_CODING
              #$con->consequence_type seems to be incorrect for sites encoding 'U'
                        if ( $con->{type} eq 'STOP_LOST' ) {
                            $con->{type} = 'NON_SYNONYMOUS_CODING';
                        }
                    }
                    else {
                        warn(      "The alignment '"
                                . $seq1->seq()
                                . "' and reference '$con->{amino_acid_reference}' amino acids do not match for sequence $con->{protein_id}."
                        );
                    }
                }
                }
                push(
                    @{ $con->{aligned_homolog_residues} },
                    $seq2->seq()
                );
		#my $test = space_to_underscore( $taxon2->short_name );
                push(
                    @{ $con->{homolog_species} },
                    space_to_underscore( $taxon2->ensembl_alias )
                );
		my $specy=space_to_underscore( $taxon2->ensembl_alias ) ;
		my $res= $seq2->seq();
		#my %results;
		if (exists $results{$specy} || (! defined $specy) ){
			next;
		}
		else {
			$results{$specy}=$res;
		}
	
		 #push( @{ $results{species} }, space_to_underscore( $taxon2->ensembl_alias ) );
		  #push( @{ $results{residue} }, $seq2->seq() );
		#print "happy";
            }
            else {
                die("Unexpected sequence ordering in alignment.");
	    }
	    
	    
	#    push (
	#	@{ $con->{orthology} },
	#	\%results
	#	);
	
	    
    } return \%results;
    

    #sort so that reference sequence is first
    #$results{species} = sort_seqs_so_reference_first( $consequence->{protein_id},
    #    $results{seqs} );
    #
    #return \%results;
    
	#my $unique_species= get_unique( \%results);
    
    
    #print "hi";
}


sub sort_seqs_so_reference_first {
    my $reference_id  = shift;
    my $array_of_seqs = shift;

    my @sorted = map { $_->[0] }

        sort { $b->[1] <=> $a->[1] }
        map {
        [   $_,
            get_sort_seqs_so_reference_first_value( $_->id, $reference_id )
        ]
        } @{$array_of_seqs};

    return \@sorted;

}

sub space_to_underscore {
    my $text = shift;
    $text =~ s/\s{1,}/_/g;
    return $text;
}

sub determine_alignment_score_change {
    my $matrix      = shift;
    my $con = shift;


    if ( scalar( @{ $con->{aligned_homolog_residues} } ) == 0 ) {
        return;
    }

    my $reference_score = get_average_similarity_score(
        $matrix,
        $con->{amino_acid_reference},
        $con->{aligned_homolog_residues}
    );
    my $variant_score = get_average_similarity_score(
        $matrix,
        $con->{amino_acid_variant},
        $con->{aligned_homolog_residues}
    );
    if (defined $variant_score && defined $reference_score) {
	$con->{alignment_score_change}
        = sprintf( "%.1f", $variant_score - $reference_score );
    }
    else {
	$con->{alignment_score_change} = undef;
    }
}
sub get_unique {
    my $list = shift;
    my %seen = ();
    my @uniq = grep { !$seen{$_}++ } @{$list};
    return \@uniq;
}

sub determine_reference_amino_acid_conservation {
    my $con = shift;

    if ( scalar( @{ $con->{aligned_homolog_residues} } ) == 0 ) {
        return;
    }

    my $reference_amino_acid_conservation = get_percent_identity_column(
        $con->{amino_acid_reference},
        $con->{aligned_homolog_residues}
    );

    if ( defined($reference_amino_acid_conservation) ) {
        $con->{reference_amino_acid_conservation}
            = sprintf( "%.1f", $reference_amino_acid_conservation );
    }
}

sub determine_reference_amino_acid_score {
    my $matrix      = shift;
    my $con = shift;

    if ( scalar( @{ $con->{aligned_homolog_residues} } ) == 0 ) {
        return;
    }

    my $reference_amino_acid_score = get_score_column(
        $matrix,
        $con->{amino_acid_reference},
        $con->{aligned_homolog_residues}
    );

    if ( defined($reference_amino_acid_score) ) {
        $con->{reference_amino_acid_score}
            = sprintf( "%.1f", $reference_amino_acid_score );
    }
}

#sub determine_context_alignments {
#    my $config     = shift;
#    my $con = shift;
#
#    if ( !( $con->{changes_protein} ) ) {
#        return;
#    }
#
#    my $context_start = $con->{altered_aa_start}
#        - 5;
#    my $context_end = $con->{altered_aa_start}
#        + 5;
#
#    if ( $context_start < 1 ) {
#        $context_start = 1;
#    }
#    if ( $context_end > length( $con->{protein_sequence} ) ) {
#        $context_end = length( $con->{protein_sequence} );
#    }
#
#    #$options->{ma} is a Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor
#    #$member is a Bio::EnsEMBL::Compara::Member
#    my $member = $config->{ma}
#        ->fetch_by_source_stable_id( 'ENSEMBLGENE', $con->{gene_id} );
#
#    #Rarely the $member object may be undef
#    if ( !defined($member) ) {
#        return;
#    }
#
#    #$options->{ha} is a Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor
#    #$homologies is a list of Bio::EnsEMBL::Compara::Homology objects
#
#    my $homologies = [];
#    
#    #obtain orthologues from all species
#    push(
#	 @{$homologies},
#	 @{  $config->{ha}->fetch_all_by_Member_method_link_type( $member,
#                    'ENSEMBL_ORTHOLOGUES' )
#	  }
#        );
#
#
#    foreach my $homology ( @{$homologies} ) {
#
#     #$homologues is an array ref of (2) Bio::EnsEMBL::Compara::Member objects
#     #$aln is Bio::SimpleAlign object
#
#        #added eval block on 2009-10-22 because
#        #$aln = $homology->get_SimpleAlign();
#        #throws exception for ENSBTAT00000060548
#        my $homologues = undef;
#        my $taxon1     = undef;
#        my $taxon2     = undef;
#        my $aln        = undef;
#
#        eval {
#            $homologues = $homology->gene_list();
#            $taxon1     = $$homologues[0]->taxon;
#            $taxon2     = $$homologues[1]->taxon;
#            $aln        = $homology->get_SimpleAlign();
#        };
#        if ($@) {
#            next;
#        }
#
#        if ( !( defined($aln) ) ) {
#            next;
#        }
#
#        #confirm that reference protein is in the alignment
#        if (!(  scalar(
#                    $aln->each_seq_with_id( $con->{protein_id} )
#                ) == 1
#            )
#            )
#        {
#            next;
#        }
#
#    #get an alignment containing the SNP-affected column and flanking sequence
#        my $sub_align_context = undef;
#
#        eval {
#            my $col_context_start
#                = $aln->column_from_residue_number(
#                $con->{protein_id},
#                $context_start );
#
#            my $col_context_end = $aln->column_from_residue_number(
#                $con->{protein_id}, $context_end );
#
#            $sub_align_context
#                = $aln->slice( $col_context_start, $col_context_end );
#        };
#        if ($@) {
#            next;
#        }
#
##confirm that there are two sequences in the subalignment
##(sequences are excluded if the alignment slice contains no residues from the sequence)
#        my $context_string = undef;
#        if (   ( defined($sub_align_context) )
#            && ( scalar( $sub_align_context->each_seq() ) == 2 ) )
#        {
#            my $seq1_context = $sub_align_context->get_seq_by_pos(1);
#            my $seq2_context = $sub_align_context->get_seq_by_pos(2);
#
#            my %seq_hash = ( reference => undef, homolog => undef );
#
#           #in all the examples I've examined, $seq1 is the reference sequence
#            if ( $seq1_context->display_id eq $con->{protein_id} ) {
#
#                $seq_hash{reference} = $seq1_context->seq();
#                $seq_hash{homolog}   = $seq2_context->seq();
#
#                push( @{ $con->{context_alignments} }, \%seq_hash );
#                push(
#                    @{ $con->{homolog_species_context} },
#                    space_to_underscore( $taxon2->short_name )
#                );
#            }
#            else {
#                die("Unexpected sequence ordering in alignment.");
#            }
#        }
#    }
#}

sub determine_context_average_percent_identity {
    my $config     = shift;
    my $con = shift;

    if (   ( !defined( $con->{context_alignments} ) )
        || ( scalar( @{ $con->{context_alignments} } ) == 0 ) )
    {
        return;
    }

    my $identity_sum   = 0;
    my $identity_count = 0;
    foreach
        my $subsequence_alignment ( @{ $con->{context_alignments} } )
    {
        my $percent_identity = get_percent_identity_between_two_seqs(
            $subsequence_alignment->{reference},
            $subsequence_alignment->{homolog}
        );
        if ( defined($percent_identity) ) {
            $identity_sum = $identity_sum + $percent_identity;
            $identity_count++;
        }
    }

    if ( $identity_count > 0 ) {
        $con->{context_conservation}
            = sprintf( "%.1f", ( $identity_sum / $identity_count ) );
    }
    return $con;
}

#compares $residue to each amino acid in $array_of_aligned
#and determines average similarity score
sub get_average_similarity_score {
    my $matrix           = shift;
    my $residue          = shift;
    my $array_of_aligned = shift;

    my $score_sum   = 0;
    my $score_count = 0;

    foreach my $aligned ( @{$array_of_aligned} ) {
        my $score = $matrix->get_entry( uc($residue), uc($aligned) );
        if ( defined($score) ) {
            $score_sum = $score_sum + $score;
            $score_count++;
        }
    }

    my $score = undef;
    if ( $score_count == 0 ) {
        $score = undef;
    }
    else {
        $score = $score_sum / $score_count;
        $score = sprintf( "%.2f", $score );
    }
    return $score;
}

#compares $residue to each amino acid in $array_of_aligned
#and determines percentage of residues in $array_of_aligned
#that are identical.
#$array_of_aligned is expected to only contain residues (no gaps)
#nonetheless gaps are skipped by this function.
#
#G vs GA will give 50
#G vs -G will give 100
#
#Is equivalent to "Cident" formula in "Correlated Mutations: A Hallmark of Phenotypic Amino Acid Substitutions"
#by Kowarsch et al., 2010

sub get_percent_identity_column {
    my $residue          = shift;
    my $array_of_aligned = shift;

    my $match_count   = 0;
    my $checked_count = 0;

    foreach my $aligned ( @{$array_of_aligned} ) {
        if ( ( is_gap($residue) ) || ( is_gap($aligned) ) ) {
            next;
        }
        $checked_count++;
        if ( uc($residue) eq uc($aligned) ) {
            $match_count++;
        }
    }
    if ( $checked_count != 0 ) {
        my $percent_identity
            = sprintf( "%.1f", ( $match_count / $checked_count ) * 100 );

        #   print "residue: $residue\n";
        #   print "aligned: " . join(',', @{$array_of_aligned}) . "\n";
        #   print "percent_identity: $percent_identity\n";

        return $percent_identity;
    }
    return undef;

}

#compares $residue to each amino acid in $array_of_aligned
#and determines sum of scores using scoring matrix.
#this value is then divided by the score obtained if all
#residues in $array_of_aligned were to match $residue
#$array_of_aligned is expected to only contain residues (no gaps)
#nonetheless gaps are skipped by this function.
#
#G vs GA will give score(G,G) + score(G,A) / 2 * score(G,G)
#G vs -G will give score(G,G) / score(G,G)
#
#Is equivalent to "Cblosum" formula in "Correlated Mutations: A Hallmark of Phenotypic Amino Acid Substitutions"
#by Kowarsch et al., 2010
sub get_score_column {
    my $matrix           = shift;
    my $residue          = shift;
    my $array_of_aligned = shift;

    my $score_sum     = 0;
    my $max_score_sum = 0;

    my $max_score = $matrix->get_entry( uc($residue), uc($residue) );
    if ( !defined($max_score) ) {
        return undef;
    }

    foreach my $aligned ( @{$array_of_aligned} ) {
        if ( ( is_gap($residue) ) || ( is_gap($aligned) ) ) {
            next;
        }

        my $score = $matrix->get_entry( uc($residue), uc($aligned) );
        if ( defined($score) ) {
            $score_sum = $score_sum + $score;
        }

        $max_score_sum = $max_score_sum + $max_score;

    }
    if ( $max_score_sum != 0 ) {
        my $score = sprintf( "%.1f", $score_sum / $max_score_sum );

        #   print "residue: $residue\n";
        #   print "aligned: " . join(',', @{$array_of_aligned}) . "\n";
        #   print "column score: $score\n";

        return $score;
    }
    return undef;

}

#determines percent identity between two sequences
#by looking for character matches.
#
#Aligning gaps are ignored
#
#--G
#--G will give 100 match
#
#-GA
#-GG will give 50 match
sub get_percent_identity_between_two_seqs {

    my $seq1 = shift;
    my $seq2 = shift;

    my @seq1 = split( //, $seq1 );
    my @seq2 = split( //, $seq2 );

    my $match_count   = 0;
    my $checked_count = 0;
    for ( my $i = 0; $i < scalar(@seq1); $i++ ) {
        if ( ( is_gap( $seq1[$i] ) ) || ( is_gap( $seq2[$i] ) ) ) {
            next;
        }

        $checked_count++;
        if ( uc( $seq1[$i] ) eq uc( $seq2[$i] ) ) {
            $match_count++;
        }
    }
    if ( $checked_count != 0 ) {
        my $percent_identity
            = sprintf( "%.1f", ( $match_count / $checked_count ) * 100 );

        #   print "seq1: $seq1\n";
        #   print "seq2: $seq2\n";
        #   print "percent_identity: $percent_identity\n";

        return $percent_identity;
    }
    return undef;
}

sub is_gap {
    my $residue = shift;
    if ( $residue =~ m/[\.\-]/ ) {
        return 1;
    }
    return 0;
}


sub determine_other_annotations {
	
   my $config = shift;
   my $con =shift;
	#Other annotations can contain a variety of key-value pairs
    my @other_annotations = ();
    
     if (( defined( $con->{entrez_gene_model_summary} ) )
        && (scalar( @{ $con->{entrez_gene_model_summary} } )
            > 0 )
        )
    {
        s/[\|\;]//g
            for ( @{ $con->{entrez_gene_model_summary} } );
        push(
            @other_annotations,
            'Gene_summary='
                . join( '|',
                @{ $con->{entrez_gene_model_summary} } )
        );
    }
    if (( defined( $con->{entrez_gene_model_gene_rif} ) )
        && (scalar( @{ $con->{entrez_gene_model_gene_rif} } )
            > 0 )
        )
    {
        s/[\|\;]//g
            for ( @{ $con->{entrez_gene_model_gene_rif} } );
        push(
            @other_annotations,
            'Function='
                . join( '|',
                @{ $con->{entrez_gene_model_gene_rif} } )
        );
    }
    
    
    my $other_an;
    if ( scalar(@other_annotations) > 0 ) {
        s/\;//g for ( @{ $con->{entrez_gene_model_phenotypes} } );
        $con->{Other_Annotations} = join( ';', @other_annotations );
	 $other_an=join(';', @other_annotations);
    }
    return $other_an;
    ########################
}

sub determine_model_annotations {
	
   my $config = shift;
   my $con =shift;
	#Model annotations can contain a variety of key-value pairs
    my @model_annotations = ();

#Phenotypic information associated with known variation at site (or orthologous site) in the model species
#using information from Ensembl
    if (   ( defined( $con->{model_phenotype_annotations} ) )
        && ( scalar( @{ $con->{model_phenotype_annotations} } ) > 0 )
        )
    {
        s/[\|\;]//g for ( @{ $con->{model_phenotype_annotations} } );
        push( @model_annotations,'Phenotypes_Position='. join( '|',@{ $con->{model_phenotype_annotations} } ) );
    }

#Phenotypes involving any region of gene from Entrez Gene entry for model gene or model orthologue
    if (   ( defined( $con->{entrez_gene_model_phenotypes} ) )
        && ( scalar( @{ $con->{entrez_gene_model_phenotypes} } ) > 0 )
        )
    {
        s/[\|\;]//g for ( @{ $con->{entrez_gene_model_phenotypes} } );
        push(
            @model_annotations,
            'Phenotypes_Gene='
                . join(
                '|', @{ $con->{entrez_gene_model_phenotypes} }
                )
        );
    }

#The number of interactions involving proteins from Entrez Gene entry for model gene or model orthologue
    if ( defined( $con->{entrez_gene_model_interactions_count} ) ) {
        push( @model_annotations,
            'Interactions_Count='
                . $con->{entrez_gene_model_interactions_count} );
    }

    #KEGG pathways from Entrez Gene entry for model gene or model orthologue
    if (( defined( $con->{entrez_gene_model_kegg_pathways} ) )
        && (scalar( @{ $con->{entrez_gene_model_kegg_pathways} } )
            > 0 )
        )
    {
        s/[\|\;]//g
            for ( @{ $con->{entrez_gene_model_kegg_pathways} } );
        push(
            @model_annotations,
            'KEGG='
                . join( '|',
                @{ $con->{entrez_gene_model_kegg_pathways} } )
        );
    }
    
    # if (( defined( $con->{entrez_gene_model_gene_rif} ) )
    #    && (scalar( @{ $con->{entrez_gene_model_gene_rif} } )
    #        > 0 )
    #    )
    #{
    #    s/[\|\;]//g
    #        for ( @{ $con->{entrez_gene_model_gene_rif} } );
    #    push(
    #        @model_annotations,
    #        'Function='
    #            . join( '|',
    #            @{ $con->{entrez_gene_model_gene_rif} } )
    #    );
    #}
    #
    ## if (( defined( $con->{entrez_gene_model_other_source_anchor} ) )
    ##    && (scalar( @{ $con->{entrez_gene_model_other_source_anchor} } )
    ##        > 0 )
    ##    )
    ##{
    ##    s/[\|\;]//g
    ##        for ( @{ $con->{entrez_gene_model_other_source_anchor} } );
    ##    push(
    ##        @model_annotations,
    ##        'Function='
    ##            . join( '|',
    ##            @{ $con->{entrez_gene_model_other_source_anchor} } )
    ##    );
    ##}
    
     if (( defined( $con->{entrez_gene_model_reactome} ) )
        && (scalar( @{ $con->{entrez_gene_model_reactome} } )
            > 0 )
        )
    {
        s/[\|\;]//g
            for ( @{ $con->{entrez_gene_model_reactome} } );
        push(
            @model_annotations,
            'Reactome_Event='
                . join( '|',
                @{ $con->{entrez_gene_model_reactome} } )
        );
    }
    #
    # if (( defined( $con->{entrez_gene_model_summary} ) )
    #    && (scalar( @{ $con->{entrez_gene_model_summary} } )
    #        > 0 )
    #    )
    #{
    #    s/[\|\;]//g
    #        for ( @{ $con->{entrez_gene_model_summary} } );
    #    push(
    #        @model_annotations,
    #        'Gene_summary='
    #            . join( '|',
    #            @{ $con->{entrez_gene_model_summary} } )
    #    );
    #}

  #Protein features from orthologues that overlap with the position of the SNP
    if ((   defined(
                $con
                    ->{uniprot_overlapping_protein_features_orthologues}
            )
        )
        && (scalar(
                @{  $con
                        ->{uniprot_overlapping_protein_features_orthologues}
                    }
            ) > 0
        )
        )
    {
        s/[\|\;]//g
            for (
            @{  $con
                    ->{uniprot_overlapping_protein_features_orthologues}
            }
            );
        push(
            @model_annotations,
            'Overlapping_Protein_Features='
                . join(
                '|',
                @{  $con
                        ->{uniprot_overlapping_protein_features_orthologues}
                    }
                )
        );
    }
    my $model_an;
    if ( scalar(@model_annotations) > 0 ) {
        s/\;//g for ( @{ $con->{entrez_gene_model_phenotypes} } );
        $con->{Model_Annotations} = join( ';', @model_annotations );
	 $model_an=join(';', @model_annotations);
    }
    return $model_an;
    ########################
}



#adds information about known model phenotypes to consequence
sub get_model_phenotypes {
    my $config     = shift;
    my $con = shift;
    my $new_vf =shift;

    my @aligned_phenotypes = ();

    if (   ( !defined( $config->{translation_adaptor} ) )
        || ( !defined( $config->{ta} ) )
        || ( !defined( $config->{va} ) )
        || ( !defined( $config->{vfa} ) )
        || ( !defined( $config->{sa} ) )
        || ( !defined( $config->{vaa} ) ) )
    {
        return;
    }

    #if input SNPs are from model then look for known variants
    #and phenotypes at that position
    if ( lc( $config->{species} ) eq lc( $config->{model} ) ) {
        my $strand=1;

        my %genomic_coords = (
            chr    => $new_vf->seq_region_name,
            start  => $new_vf->start,
            end    => $new_vf->start,
            strand => $strand
        );
        push(
            @aligned_phenotypes,
            @{  get_model_phenotypes_from_model_genomic_coords( $config,
                    \%genomic_coords )
                }
        );

    }

    #otherwise transform non-model coordinate to model
    else {
        if (   ( !defined( $con->{protein_id} ) )
            || ( !defined( $con->{altered_aa_start} ) )
            || ( !defined( $con->{altered_aa_end} ) ) )
        {
            return;
        }

        foreach my $model_orthologous_residue (
            @{ $con->{aligned_model_orthologue_proteins} } )
        {
            my $genomic_coords= get_model_genomic_coords_of_model_amino_acid( $config, $model_orthologous_residue );
            if ( defined($genomic_coords) ) {
                push(
                    @aligned_phenotypes,
                    @{  get_model_phenotypes_from_model_genomic_coords($config, $genomic_coords )
                        }
                );
            }
        }
    }

    my $unique_phenotypes = get_unique( \@aligned_phenotypes );
    push(
        @{ $con->{model_phenotype_annotations} },
        @{$unique_phenotypes}
    );
}


#returns array of phenotypes associated with model genomic coordinates
sub get_model_phenotypes_from_model_genomic_coords {
    my $config        = shift;
    my $genomic_coords = shift;

    my @phenotypes = ();

    my $slice = $config->{sa}->fetch_by_region(
        undef,
        $genomic_coords->{chr},
        $genomic_coords->{start},
        $genomic_coords->{end},
        $genomic_coords->{strand}
    );

    my $vfs = $config->{vfa}
        ->fetch_all_by_Slice($slice);

    foreach my $vf ( @{$vfs} ) {

        my $variation_name = $vf->variation_name();
        my $variation      = $config->{va}->fetch_by_name($variation_name);
        my $variation_class = $vf->var_class();

        my $vas = $config->{vaa}->fetch_all_by_Variation($variation);

        foreach my $va ( @{$vas} ) {

#           print $va->phenotype_name(), $va->phenotype_description(), $va->source_name(), $va->study_type(), $va->local_stable_id(),"\n";
            my @phenotype_descriptors = ();
            if ( defined( $va->source_name() ) ) {
                push( @phenotype_descriptors,
                    'Source: ' . $va->source_name() );
            }
            if ( defined( $va->phenotype_name() ) ) {
                push( @phenotype_descriptors,
                    'Phenotype name: ' . $va->phenotype_name() );
            }
            if ( defined( $va->phenotype_description() ) ) {
                push( @phenotype_descriptors,
                    'Description: ' . $va->phenotype_description() );
            }
            if ( defined( $va->variation_names() ) ) {
                push( @phenotype_descriptors,
                    'Variation names: ' . $va->variation_names() );
            }
            if ( scalar(@phenotype_descriptors) > 0 ) {
                push( @phenotypes, join( ' ', @phenotype_descriptors ) );
            }
        }
    }
    return \@phenotypes;
}


#returns hash containing chromosome, start, end, and strand of genomic sequence encoding
#amino acid described by $model_translation_id, $amino_acid_start, $amino_acid_end
sub get_model_genomic_coords_of_model_amino_acid {
    my $config                   = shift;
    my $model_orthologous_residue = shift;

    my $model_translation_id = $model_orthologous_residue->{id};
    my $amino_acid_start     = $model_orthologous_residue->{start};
    my $amino_acid_end       = $model_orthologous_residue->{end};

    #obtain the model protein translation object
    my $model_translation = $config->{translation_adaptor}
        ->fetch_by_stable_id($model_translation_id);
    my $model_transcript = $config->{ta}
        ->fetch_by_translation_stable_id($model_translation_id);

    #    print "\$amino_acid_start is $amino_acid_start\n";
    #    print "\$amino_acid_end is $amino_acid_end\n";

    #determine the position on the transcript
    my $trmapper = $model_transcript->get_TranscriptMapper();

    my $model_slice = $config->{sa}->fetch_by_transcript_stable_id( $model_transcript->stable_id );

    my $model_chromosome_name = $model_slice->seq_region_name();
    my $model_slice_start     = $model_slice->start();
    my $model_slice_end       = $model_slice->end();
    my $model_slice_strand    = $model_slice->strand();

#a list of Bio::EnsEMBL::Mapper::Gap and Bio::EnsEMBL::Mapper::Coordinate objects
    my @model_genomic_coords
        = $trmapper->pep2genomic( $amino_acid_start, $amino_acid_end );

    #not sure how to handle more complicated coordinates, so skip
    if ( scalar(@model_genomic_coords) > 1 ) {
        return;
    }

    my $model_aligned_aa_start_genomic = $model_genomic_coords[0]->start;
    my $model_aligned_aa_end_genomic   = $model_genomic_coords[0]->end;
    my $model_aligned_strand_genomic   = $model_genomic_coords[0]->strand;

#check that the slice start and end is consistent with the position obtained using pep2genomic
    if (   ( !( $model_slice_start <= $model_aligned_aa_start_genomic ) )
        || ( !( $model_slice_end >= $model_aligned_aa_end_genomic ) ) )
    {
        return;
    }

    my %model_genomic_cords = (
        chr    => $model_chromosome_name,
        start  => $model_aligned_aa_start_genomic,
        end    => $model_aligned_aa_end_genomic,
        strand => $model_aligned_strand_genomic
    );
    return \%model_genomic_cords;

}

#Get overlapping protein features orthologous proteins, from UniProt
sub get_overlapping_protein_features_from_UniProt_orthologues {
    my $config     = shift;
    my $con = shift;
    my $parsed =shift;

    foreach my $orthologous_protein (
        @{ $con->{aligned_model_orthologue_proteins} } )
    {

        if (   ( !defined( $orthologous_protein->{uniprot_id} ) )
            || ( !defined( $orthologous_protein->{sequence} ) )
            || ( !defined( $orthologous_protein->{start} ) )
            || ( !defined( $orthologous_protein->{end} ) ) )
        {
            next;
        }

        foreach my $id ( @{ $orthologous_protein->{uniprot_id} } ) {

            my $record = get_uniprot_record(
                script => $config->{uniprot_script},
                id     => $id,
                db     => 'uniprotkb',
                format => 'SWISS',
                style  => 'raw'
            );
	    
            my $overlapping_features
                = get_overlapping_features_from_uniprot_record(
                record            => $record,
                expected_sequence => $orthologous_protein->{sequence},
                region_start      => $orthologous_protein->{start},
                region_end        => $orthologous_protein->{end},
                id                => $id
                );

            if ( defined($overlapping_features) ) {
                push(
                    @{  $con->{uniprot_overlapping_protein_features_orthologues}
                        },
                    @{$overlapping_features}
                );
            }
        }
    }
    $con->{uniprot_overlapping_protein_features_orthologues}
        = get_unique(
        $con->{uniprot_overlapping_protein_features_orthologues} );
}

#Get overlapping protein features from UniProt
sub get_overlapping_protein_features_from_UniProt {
    my $config     = shift;
    my $con = shift;
    my $parsed = shift;

    if (   
        ( !defined( $con->{protein_sequence} ) )
        || ( !defined( $con->{altered_aa_start} ) )
        || ( !defined( $con->{altered_aa_end} ) ) )
    {
        return;
    }

    foreach my $id ( @{ $con->{uniprot_id} } ) {

        my $record = get_uniprot_record(
            script => $config->{uniprot_script},
            id     => $id,
            db     => 'uniprotkb',
            format => 'SWISS',
            style  => 'raw'
        );
	
	my $entry = parse_uniprot_flatfile($record, $con);
	print $entry->{protein_length};

#my $extra;
    foreach my $key (keys %$entry) {
	if ( ( defined( $entry->{$key} ) ) && scalar( @{ $entry->{$key} } ) > 0)  {
		$con->{uniprot_parsed} .= uc($key)."=".join('|', @{ $entry->{$key} }).';';
		#$extra .= uc($key)."=".join('|', @{ $entry->{$key} }).';';
	}
    }

        my $overlapping_features
            = get_overlapping_features_from_uniprot_record(
            record            => $record,
            expected_sequence => $con->{protein_sequence},
            region_start      => $con->{altered_aa_start},
            region_end        => $con->{altered_aa_end},
            id                => $id
            );
        if ( defined($overlapping_features) ) {
            push(
                @{ $con->{uniprot_overlapping_protein_features} },
                @{$overlapping_features}
            );
        }
    }
    $con->{uniprot_overlapping_protein_features}
        = get_unique( $con->{uniprot_overlapping_protein_features} );
}

sub get_overlapping_features_from_uniprot_record {
    my %args = (@_);

    if ( !defined( $args{record} ) ) {
        return undef;
    }
  
    my $sequence;
    if ( $args{record} =~ m/^SQ\s{3}SEQUENCE[^\n]+\n([^\/]+)\/\//m ) {
        $sequence = $1;
        $sequence =~ s/\s//g;
    }
    if ( !defined($sequence) ) {

        #sometimes the accession obtained from Ensembl
        #doesn't work for retrieval, e.g. ANK36_HUMAN
        return undef;
    }

    if ( uc($sequence) ne uc( $args{expected_sequence} ) ) {

        #sometimes the sequences don't match
        #do not transfer features if this is the case
        return undef;
    }

    #obtain features
    #FT   CHAIN         1    190       Selenoprotein S.
    #FT                                /FTId=PRO_0000318650.
    #FT   TRANSMEM     28     48       Helical; (Potential).
    #FT   NON_STD     189    189       Selenocysteine (By similarity).
    #FT   MOD_RES     140    140       Phosphoserine (By similarity).

    my @features = ();
    while ( $args{record} =~ m/^FT\s{3}(\S+)\s+(\d+)\s+(\d+)\s+(\S[^\n]+)/mg )
    {
        my $key         = $1;
        my $start       = $2;
        my $end         = $3;
        my $description = $4;

        if (overlaps(
                s1 => $args{region_start},
                e1 => $args{region_end},
                s2 => $start,
                e2 => $end
            ) )
        {
            push( @features, "$key:$start:$end:$description" );
        }
    }
            if ( scalar(@features) > 0 ) {
        return \@features;
    }
    return undef;
}

#Get UniProt record from EBI
sub get_uniprot_record {
    my %args = (@_);

    #create temp file for output
    my $tmp_output          = new File::Temp();
    my $tmp_output_filename = $tmp_output->filename;

    my $command
        = 'perl '
        . $args{script}
        . " -i '$args{id}'"
        . " -o $tmp_output_filename"
        . " -d '$args{db}'"
        . " -f '$args{format}'"
        . " -s '$args{style}'";

    my $result = system($command);
    if ( $result != 0 ) {
        die("The following command failed: '$command'\n");
    }

    close($tmp_output) or die("Cannot close file : $!");

    local $/ = undef;
    open( my $FILE, '<', $tmp_output_filename )
        or die("Cannot open file '$tmp_output_filename': $!");
	
    my $output = <$FILE>;


    close($FILE) or die("Cannot close file : $!");
    return $output;
}


#Get information about model orthologue from NCBI Gene database.
#Rationale is that model genes have much more detailed annotation.
sub get_model_orthologue_info_from_NCBI {
    my $config     = shift;
    my $con = shift;
    my $parsed = shift;



    #both of the following can be used to retrieve Entrez Gene records
    my $model_entrez_gene_id;
    my $model_refseq_peptide_accession;

    if ( lc( $config->{species} ) eq lc( $config->{model} ) ) {
        if ( defined( $con->{entrez_gene_id} ) ) {
            $model_entrez_gene_id = $con->{entrez_gene_id};
        }
        else {
            return;
        }
    }
    else {
        if ( !defined( $con->{gene_id} ) ) {
            return;
        }
		
		
        #$options->{ma} is a Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor
        #$member is a Bio::EnsEMBL::Compara::Member
        my $member = $config->{ma}->fetch_by_source_stable_id( 'ENSEMBLGENE',
            $con->{gene_id} );

        #Rarely the $member object may be undef
        if ( !defined($member) ) {
            return;
        }

        #$options->{ha} is a Bio::EnsEMBL::Compara::DBSQL::HomologyAdaptor
        #$homologies is a list of Bio::EnsEMBL::Compara::Homology objects

        my $homologies = [];
        push(
            @{$homologies},
            @{  $config->{ha}->fetch_all_by_Member_paired_species( $member,
                    $config->{model}, ['ENSEMBL_ORTHOLOGUES'] )
                }
        );

        foreach my $homology ( @{$homologies} ) {

     #$homologues is an array ref of (2) Bio::EnsEMBL::Compara::Member objects

            my $homologues     = undef;
            my $homologue_gene = undef;

            eval {
                $homologues     = $homology->gene_list();
                $homologue_gene = $$homologues[1]->get_Gene;
            };
            if ($@) {
                next;
            }

            my $xrefs = $homologue_gene->get_all_DBLinks();

            #try to obtain Entrez Gene ID first
            foreach my $xref ( @{$xrefs} ) {
                if ( $xref->dbname() eq 'EntrezGene' ) {
                    $model_entrez_gene_id = $xref->primary_id();
                }
            }
            if ( defined($model_entrez_gene_id) ) {
                last;
            }

            #try to obtain RefSeq protein ID second
            foreach my $xref ( @{$xrefs} ) {
                if ( $xref->dbname() eq 'RefSeq_peptide' ) {
                    $model_refseq_peptide_accession = $xref->primary_id();
                }
            }
            if ( defined($model_refseq_peptide_accession) ) {
                last;
            }

        }

    }
    if (defined ($parsed->{$con->{entrez_gene_id} })) {
    $con->{entrez_gene_model_kegg_pathways} = $parsed->{$con->{entrez_gene_id} }->{kegg};
    $con->{entrez_gene_model_phenotypes}    = $parsed->{$con->{entrez_gene_id} }->{phenotypes};
    $con->{entrez_gene_model_interactions_count}= $parsed->{$con->{entrez_gene_id} }->{interaction_count};
    $con->{entrez_gene_model_reactome} = $parsed->{$con->{entrez_gene_id} }->{reactome};
    $con->{entrez_gene_model_summary} = $parsed->{$con->{entrez_gene_id} }->{entrez_gene_summary};
    #$con->{entrez_gene_model_other_source_anchor} = $parsed->{other_source_anchor};
    $con->{entrez_gene_model_gene_rif} = $parsed->{$con->{entrez_gene_id} }->{gene_rif};
    return;
					}
    my $xml;
    if ( defined($model_entrez_gene_id) ) {
        $xml = get_ncbi_record(
            script => $config->{ncbi_script},
            db     => 'gene',
            query  => $model_entrez_gene_id . '[UID]'
        );
    }
    
    #Use BRCA for debugging purposes
     #if ( defined($model_entrez_gene_id) ) {
     #   $xml = get_ncbi_record(
     #       script => $config->{ncbi_script},
     #       db     => 'gene',
     #       query  => 'BRCA2'
     #   );
     #}
     
       elsif ( defined($model_refseq_peptide_accession) ) {
        $xml = get_ncbi_record(
            script => $config->{ncbi_script},
            db     => 'gene',
            query  => $model_refseq_peptide_accession . '[Protein Accession]'
        );
    }
       
    $parsed->{$con->{entrez_gene_id} } = get_info_from_ncbi_gene_record($xml);
    print "parsing ".$con->{entrez_gen_id}."\n";

    $con->{entrez_gene_model_kegg_pathways} = $parsed->{$con->{entrez_gene_id} }->{kegg};
    $con->{entrez_gene_model_phenotypes}    = $parsed->{$con->{entrez_gene_id} }->{phenotypes};
    $con->{entrez_gene_model_interactions_count}= $parsed->{$con->{entrez_gene_id} }->{interaction_count};
    $con->{entrez_gene_model_reactome} = $parsed->{$con->{entrez_gene_id} }->{reactome};
    $con->{entrez_gene_model_summary} = $parsed->{$con->{entrez_gene_id} }->{entrez_gene_summary};
    #$con->{entrez_gene_model_other_source_anchor} = $parsed->{other_source_anchor};
    $con->{entrez_gene_model_gene_rif} = $parsed->{$con->{entrez_gene_id} }->{gene_rif};
    

}

sub get_info_from_ncbi_gene_record {
    my $xml = shift;
    my %results
        = ( kegg => undef, phenotypes => undef, interaction_count => undef , reactome => undef, entrez_gene_summary => undef, other_source_anchor=> undef, gene_rif => undef);

    if ( !defined($xml) ) {
        return \%results;
    }

    #KEGG pathways
    #<Gene-commentary_text>KEGG pathway: Apoptosis</Gene-commentary_text>
    my $kegg_re
        = "\Q<Gene-commentary_text>KEGG pathway: \E(.+)\Q</Gene-commentary_text>\E";

    $results{kegg} = get_unique_re_matches( $kegg_re, $xml );
    
    my $reactome_re
        ="\Q<Gene-commentary_text>Reactome Event:\E(.+)\Q</Gene-commentary_text>\E";
	
	$results{reactome}=get_unique_re_matches( $reactome_re, $xml);
	
	
    my $entrez_gene_summary_re
        ="\Q<Entrezgene_summary>\E(.+)\Q</Entrezgene_summary>\E";
	
	$results{entrez_gene_summary}=get_unique_re_matches( $entrez_gene_summary_re, $xml);


#	my $gene_rif_re    
#        = "\Q<Gene-commentary_type value=\"generif\">18</Gene-commentary_type>\E"
#        . '[\s\n]*'
#        . "\Q<Gene-commentary_text>\E(.+)\Q</Gene-commentary_text>\E";
#
#    $results{gene_rif} = get_unique_re_matches( $gene_rif_re, $xml );


    #phenotypes
    #<Gene-commentary_type value="phenotype">19</Gene-commentary_type>
    #<Gene-commentary_text>Colorectal cancer</Gene-commentary_text>
    my $phenotypes_re
        = "\Q<Gene-commentary_type value=\"phenotype\">19</Gene-commentary_type>\E"
        . '[\s\n]*'
        . "\Q<Gene-commentary_text>\E(.+)\Q</Gene-commentary_text>\E";

    $results{phenotypes} = get_unique_re_matches( $phenotypes_re, $xml );

    #interactions
    #<Gene-commentary_heading>Interactions</Gene-commentary_heading>
    #to
    #<Gene-commentary_heading>Pathways</Gene-commentary_heading>
    my $interactions_section;
    if ( $xml
        =~ m/<Gene\-commentary_heading>Interactions<\/Gene\-commentary_heading>([\s\S]+)<Gene\-commentary_heading>Pathways<\/Gene\-commentary_heading>/
        )
    {
        $interactions_section = $1;
    }
    if ( defined($interactions_section) ) {
        my $interactions_re
            = "\Q<Gene-commentary_type value=\"generif\">18</Gene-commentary_type>\E"
            . '[\s\n]*'
            . "\Q<Gene-commentary_text>\E(.+)\Q</Gene-commentary_text>\E";

        my $interactions = get_unique_re_matches( $interactions_re,
            $interactions_section );
        if ( defined($interactions) ) {
            $results{interaction_count} = scalar( @{$interactions} );
        }

    }
    return \%results;
}

sub get_ncbi_record {
    my %args = (@_);

    my $tmp_output          = new File::Temp();
    my $tmp_output_filename = $tmp_output->filename;

    my $command
        = 'perl '
        . $args{script}
        . " -q '$args{query}'"
        . " -o $tmp_output_filename"
        . " -d '$args{db}'"
        . " -r xml" . " -m 1";

    my $result = system($command);
    if ( $result != 0 ) {
        die("The following command failed: '$command'\n");
    }

    close($tmp_output) or die("Cannot close file : $!");

    local $/ = undef;
    open( my $FILE, '<', $tmp_output_filename )
        or die("Cannot open file '$tmp_output_filename': $!");

    my $output = <$FILE>;

    close($FILE) or die("Cannot close file : $!");
    return $output;
}

sub get_unique_re_matches {
    my $re      = shift;
    my $text    = shift;
    my @matches = ();
    while ( $text =~ m/$re/g ) {
        push( @matches, $1 );
    }

    if ( ( scalar(@matches) ) == 0 ) {
        return undef;
    }
    return get_unique( \@matches );
}

sub parse_uniprot_flatfile {
    my $record = shift;
    my $con = shift;
    
    local $/ = "\n//\n";
    my $fullParse=0;
    
    my %results = (
                'description'   => undef,
                'modified_positions'  => undef,
                'rna_editing' => undef,
                'development'  => undef, 
                'dis_phenotype'  =>undef,
                'domains'  => undef, 
                'enzyme_reg'  => undef, 
                'function'  => undef, 
                'induction'  => undef, 
                'interactions'  => undef, 
                'allergens'  => undef, 
                'misc'  => undef, 
                'pathway'  => undef, 
                'polymorphisms'  => undef, 
                'ptms'  => undef, 
                'similarity'  => undef,
                'keywords'   => undef,
                'tissue'  => undef, 
                'subunit'  => undef,
                'subcell_location' => undef
            );

#while (<>){
  # Read the entry
  my  $entry = SWISS::Entry->fromText($record, $fullParse);
  
  $con->{protein_length} = $entry->IDs->length;
  
  foreach my $de ($entry->DEs->elements) {
   push @{$results{description}}, $de->text;
  }
  # Print all keywords
  foreach my $kw ($entry->KWs->elements) {
    push @{$results{keywords}}, $kw->text;
  }
  my @CCs = $entry->CCs->elements();

  for my $CC (@CCs) {
    if ($CC -> topic eq 'RNA EDITING') {
        #print $CC->toString;
        push @{$results{rna_editing}}, $CC->note;
        #my $term = $CC->term;
        #my $el = $CC->list->elements;
        #print $CC->text;
        #push @modified_positions, $CC->list;
    }
    elsif ($CC -> topic eq 'DEVELOPMENTAL STAGE') {
         push @{$results{development}}, $CC->comment;
    }
    elsif ($CC -> topic eq 'DISEASE') {
        push @{$con->{disease}}, $CC->comment;   
    }
    elsif ($CC -> topic eq 'DISRUPTION PHENOTYPE') {
        push @{$results{dis_phenotype}}, $CC->comment;
    }
    elsif ($CC -> topic eq 'ALLERGENS') {
        push @{$results{allergens}}, $CC->comment;
    }
    elsif ($CC -> topic eq 'DOMAIN') {
        push @{$results{domains}}, $CC->comment;   
    }
    elsif ($CC -> topic eq 'FUNCTION') {
        push @{$results{function}}, $CC->comment;
    }
    elsif ($CC -> topic eq 'PATHWAY') {
        push @{$results{pathway}}, $CC->comment;
    }
    elsif ($CC -> topic eq 'POLYMORPHISM') {
        push @{$results{polymorphisms}}, $CC->comment;   
    }
    elsif ($CC -> topic eq 'SIMILARITY') {
        push @{$results{similarity}}, $CC->comment;
    }
     elsif ($CC -> topic eq 'TISSUE SPECIFICITY') {
        push @{$results{tissue}}, $CC->comment;   
    }
    elsif ($CC -> topic eq 'PTM') {
        push @{$results{ptms}}, $CC->comment;
    }   
    else { next }
  }
#}
 
  return \%results;
	
}

sub overlaps {
    my %args = (@_);

    if (   ( $args{s1} >= $args{s2} )
        && ( $args{e1} <= $args{e2} ) )
    {
        return 1;
    }
    elsif (( $args{s1} <= $args{s2} )
        && ( $args{e1} >= $args{e2} ) )
    {
        return 1;
    }
    elsif (
        ( $args{s1} <= $args{s2} )
        && (   ( $args{e1} <= $args{e2} )
            && ( $args{e1} >= $args{s2} ) )
        )
    {
        return 1;
    }
    elsif (
        ( $args{e1} >= $args{e2} )
        && (   ( $args{s1} >= $args{s2} )
            && ( $args{s1} <= $args{e2} ) )
        )
    {
        return 1;
    }
    return 0;
}
