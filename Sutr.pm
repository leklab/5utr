=head1 LICENSE

Implemented from CADD plugin @ https://github.com/Ensembl/VEP_plugins/blob/release/96/CADD.pm

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 implemented by sander.pajusalu@yale.edu
 For Ensembl:
 Ensembl <http://www.ensembl.org/info/about/contact/index.html>
    
=cut

=head1 NAME

 5utr/Sutr

=head1 SYNOPSIS

./vep -i variants.vcf.gz --cache --offline --assembly GRCh37 --dir_plugin /path/to/Sutr_dir --no_stats --plugin Sutr,deltas_vep.tsv.gz -o out.file.txt

=head1 DESCRIPTION

 A VEP plugin that retrieves 5'utr annotation scores for variants from one or more
 tabix-indexed 5utr data files.
 
 Please cite the 5utr and data source publications alongside the VEP if you use this resource:
 
 
 The tabix utility must be installed in your path to use this plugin. The Sutr
 data files can be downloaded from github repo.
 
=cut

package Sutr;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);

  $self->get_user_params();

  return $self;
}

sub feature_types {
  return ['Feature','Intergenic'];
}

sub get_header_info {
  my $self = shift;
  return {
    delta_TE => 'predicted Transcriptional Efficiency change from reference',
    delta_rG4 => 'predicted G4 MFE change from reference',
    delta_dsRNA => 'predicted dsRNA MFEchange from reference'
  }
}

sub run {
  my ($self, $tva) = @_;
  
  my $vf = $tva->variation_feature;
  
  # get allele, reverse comp if needed
  my $allele = $tva->variation_feature_seq;
  reverse_comp(\$allele) if $vf->{strand} < 0;
  
  return {} unless $allele =~ /^[ACGT-]+$/;

  my @data =  @{$self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})};

  foreach (@data) {
    my $matches = get_matched_variant_alleles(
      {
        ref    => $vf->ref_allele_string,
        alts   => $vf->alt_alleles,
        pos    => $vf->{start},
        strand => $vf->strand
      },
      {
       ref  => $_->{ref},
       alts => [$_->{alt}],
       pos  => $_->{start},
      }
    );
    return $_->{result} if (@$matches);
  }
  return {};
}

sub parse_data {
  my ($self, $line) = @_;
  my ($c, $s, $ref, $alt, $dTE, $dG4, $dds) = split /\t/, $line;

  # do VCF-like coord adjustment for mismatched subs
  my $e = ($s + length($ref)) - 1;
  if(length($alt) != length($ref)) {
    my $first_ref = substr($ref, 0, 1);
    my $first_alt = substr($alt, 0, 1);
    if ($first_ref eq $first_alt) {
      $s++;
      $ref = substr($ref, 1);
      $alt = substr($alt, 1);
      $ref ||= '-';
      $alt ||= '-';
    }
  }
  return {
    ref => $ref,
    alt => $alt,
    start => $s,
    end => $e,
    result => {
      delta_TE   => $dTE,
      delta_rG4 => $dG4,
      delta_dsRNA => $dds
    }
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
