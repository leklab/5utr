=head1 LICENSE

This script was derived from dbSNFP.pm script https://raw.githubusercontent.com/Ensembl/VEP_plugins/release/97/dbNSFP.pm version: c3b1143 on 24 Apr

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

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>
    
=cut

=head1 NAME

 Sutr

=head1 SYNOPSIS

 mv Sutr.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin Sutr,/path/to/Sutr.gz,col1,col2
 ./vep -i variations.vcf --plugin Sutr,'consequence=ALL',/path/to/Sutr.gz,col1,col2
 ./vep -i variations.vcf --plugin Sutr,'consequence=3_prime_UTR_variant&intron_variant',/path/to/Sutr.gz,col1,col2

=head1 DESCRIPTION

 A VEP plugin that retrieves data for 5'UTR variants from a tabix-indexed
 Sutr file.
 
 You must have the Bio::DB::HTS module or the tabix utility must be installed
 in your path to use this plugin.

 When running the plugin you must list at least one column to retrieve from the
 Sutr file, specified as parameters to the plugin e.g.
 
 --plugin Sutr,/path/to/Sutr.gz,LRT_score,GERP++_RS

 You may include all columns with ALL; this fetches a large amount of data per
 variant!:
 
 --plugin Sutr,/path/to/Sutr.gz,ALL
 
 Tabix also allows the data file to be hosted on a remote server. This plugin is
 fully compatible with such a setup - simply use the URL of the remote file:
 
 --plugin Sutr,http://my.files.com/Sutr.gz,col1,col2

 Replacement logic was remained in the script, also not relevant with the current version.
 The plugin replaces occurrences of ';' with ',' and '|' with '&'. However, some
 data field columns, e.g. Interpro_domain, use the replacement characters. We
 added a file with replacement logic for customising the required replacement
 of ';' and '|' in dbNSFP data columns. In addition to the default replacements
 (; to , and | to &) users can add customised replacements. Users can either modify
 the file dbNSFP_replacement_logic in the VEP_plugins directory or provide their own
 file as second argument when calling the plugin:

 --plugin Sutr,/path/to/Sutr.gz,/path/to/dbNSFP_replacement_logic,LRT_score,GERP++_RS
 
 If the Sutr README file is found in the same directory as the data file,
 column descriptions will be read from this and incorporated into the VEP output
 file header.
 
=cut

package Sutr;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

my %INCLUDE_SO = map {$_ => 1} qw(missense_variant stop_lost stop_gained start_lost 5_prime_UTR_variant);

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  $self->expand_left(0);
  $self->expand_right(0);
  my $index = 0;
  $self->{consequence} = 'filter';
  if ($self->params->[$index] =~ /^consequence=/) {
    # parse consequences
    my $consequences = $self->params->[$index];  
    $consequences =~ s/consequence=//;
    if (uc $consequences eq 'ALL') {
      $self->{consequence} = 'ALL';
    } else {
      %INCLUDE_SO = map {$_ => 1} split/&/, $consequences;
    }
    $index++;
  }
  
  # get Sutr file
  my $file = $self->params->[$index];

  $self->add_file($file);
  
  # get headers
  open HEAD, "tabix -fh $file 1:1-1 2>&1 | ";
  while(<HEAD>) {
    next unless /^\#/;
    chomp;
    $_ =~ s/^\#//;
    $self->{headers} = [split];
  }
  close HEAD;
  
  die "ERROR: Could not read headers from $file\n" unless defined($self->{headers}) && scalar @{$self->{headers}};
  
  # check ALT and transcript headers
  foreach my $h(qw(ALT transcript)) {
    die "ERROR: Could not find required column $h in $file\n" unless grep {$_ eq $h} @{$self->{headers}};
  }

  # check if 2nd argument is a file that specifies replacement logic
  # read replacement logic 
  $index++;
  my $replacement_file = $self->params->[$index];
  if (defined $replacement_file && -e $replacement_file) {
    $self->add_replacement_logic($replacement_file);  
    $index++;
  } else {
    $self->add_replacement_logic();  
  } 
 
  # get required columns
  while(defined($self->params->[$index])) {
    my $col = $self->params->[$index];
    if($col eq 'ALL') {
      $self->{cols} = {map {$_ => 1} @{$self->{headers}}};
      last;
    }
    die "ERROR: Column $col not found in header for file $file. Available columns are:\n".join(",", @{$self->{headers}})."\n" unless grep {$_ eq $col} @{$self->{headers}};
    
    $self->{cols}->{$self->params->[$index]} = 1;
    $index++;
  }
  
  die "ERROR: No columns selected to fetch. Available columns are:\n".join(",", @{$self->{headers}})."\n" unless defined($self->{cols}) && scalar keys %{$self->{cols}};
  
  return $self;
}

sub feature_types {
  return ['Transcript'];
}

sub get_header_info {
  my $self = shift;

  if(!exists($self->{_header_info})) {

    # look for readme
    my $file_dir = $self->files->[0];

    my %rm_descs;

    # won't work for remote
    if($file_dir !~ /tp\:\/\//) {

      # get just dir
      $file_dir =~ s/\/[^\/]+$/\//;

      if(opendir DIR, $file_dir) {
        my ($readme_file) = grep {/sutr\.readme\.txt/i} readdir DIR;
        closedir DIR;

        if(open RM, $file_dir.$readme_file) {
          my ($col, $reading);

          # parse Sutr readme
          # relevant lines look like:
          #
          # 1   column1_name: description blah blah
          #     blah blah blah
          # 2   column2_name: description blah blah
          #     blah blah blah

          while(<RM>) {
            chomp;
            s/\r$//g;

            if(/^\d+\s/) {
              $reading = 1;

              m/^\d+\s+(.+?)\:\s+(.+)/;
              $col = $1;

              $rm_descs{$col} = '(from Sutr) '.$2 if $col && $2;
            }
            elsif($reading && /\w/) {
              s/^\s+//;
              $rm_descs{$col} .= ' '.$_;
            }
            else {
              $reading = 0;
            }
          }

          close RM;

          # remove multiple spaces
          $rm_descs{$_} =~ s/\s+/ /g for keys %rm_descs;
        }
      }
    }

    $self->{_header_info} = {map {$_ => $rm_descs{$_} || ($_.' from Sutr file')} keys %{$self->{cols}}};
  }
  
  return $self->{_header_info};
}

sub run {
  my ($self, $tva) = @_;
  
  # only for missense variants
  if ($self->{consequence} eq 'filter') {
    return {} unless grep {$INCLUDE_SO{$_->SO_term}} @{$tva->get_all_OverlapConsequences};
  }
  
  my $vf = $tva->variation_feature;
  
  return {} unless $vf->{start} eq $vf->{end};
  
  # get allele, reverse comp if needed
  my $allele = $tva->variation_feature_seq;
  reverse_comp(\$allele) if $vf->{strand} < 0;
  
  return {} unless $allele =~ /^[ACGT]$/;
  
  # get transcript stable ID
  my $tr_id = $tva->transcript->stable_id;

  my $data;
  my $pos;

  my $assembly = $self->{config}->{assembly};
  my $chr = ($vf->{chr} =~ /MT/i) ? 'M' : $vf->{chr};
  foreach my $tmp_data(@{$self->get_data($chr, $vf->{start} - 1, $vf->{end})}) {
    # compare allele and transcript
    if ($assembly eq 'GRCh37') {
      if (exists $tmp_data->{'pos'}) {
        $pos = $tmp_data->{'pos'}
      } else {
        die "Sutr file does not contain required columns (pos) to use with GRCh37";
        }
    }
    next unless
      $pos == $vf->{start} &&
      defined($tmp_data->{ALT}) &&
      $tmp_data->{ALT} eq $allele;
    # make a clean copy as we're going to edit it
    %$data = %$tmp_data;

    #convert data with multiple transcript values
     if($data->{transcript} =~ m/\,/) {

       # find the "index" of this transcript
       my @tr_ids = split(',', $data->{transcript});
       my $tr_index;

       for my $i(0..$#tr_ids) {
         $tr_index = $i;
         last if $tr_ids[$tr_index] =~ /^$tr_id(\.\d+)?$/;
       }

       next unless defined($tr_index);

       # now alter other fields
       foreach my $key(keys %$data) {
         if($data->{$key} =~ m/\,/) {
           my @split = split(',', $data->{$key});
           #die("ERROR: Transcript index out of range") if $tr_index > $#split;
           $data->{$key} = $split[$tr_index];
         } 
       }
     }
    last;
  }
  
  return {} unless scalar keys %$data;
  
  # get required data
  my @from = @{$self->{replacement}->{default}->{from}};
  my @to = @{$self->{replacement}->{default}->{to}};

  my %return;
  foreach my $colname (keys %$data) {
    next if(!defined($self->{cols}->{$colname}));
    next if($data->{$colname} eq '.');

    my @from = @{$self->{replacement}->{default}->{from}};
    my @to   = @{$self->{replacement}->{default}->{to}};
    @from    = @{$self->{replacement}->{$colname}->{from}} if (defined $self->{replacement}->{$colname});
    @to      = @{$self->{replacement}->{$colname}->{to}} if (defined $self->{replacement}->{$colname});
    for my $i (0 .. $#from) {
      $data->{$colname} =~ s/\Q$from[$i]\E/$to[$i]/g;
    }
    $return{$colname} = $data->{$colname};
  }
  
  return \%return;
}

sub parse_data {
  my ($self, $line) = @_;

  $line =~ s/\r$//g;

  my @split = split /\t/, $line;
  
  # parse data into hash of col names and values
  my %data = map {$self->{headers}->[$_] => $split[$_]} (0..(scalar @{$self->{headers}} - 1));

  return \%data;
}

sub get_start {  
  return $_[1]->{'pos(1-based)'};
}

sub get_end {
  return $_[1]->{'pos(1-based)'};
}

sub add_replacement_logic {
  my $self = shift;
  my $file = shift;
  $file ||= 'replacement_logic';
  if (! -e $file) {
    $self->{replacement}->{default}->{from} = [';', '|'];
    $self->{replacement}->{default}->{to} = [',', '&'];
  } else {
    open FILE, $file;
    while(<FILE>) {
      chomp;
      next if /^colname/;
      my ($colname, $from, $to) = split/\s+/;
      die ("ERROR: 3 values separated by whitespace are required: colname from to.") if(!($colname && $from && $to));
      push @{$self->{replacement}->{$colname}->{from}}, $from;
      push @{$self->{replacement}->{$colname}->{to}}, $to;
    }
    close FILE;
    die("ERROR: No default replacement logic has been specified.\n") if (!defined $self->{replacement}->{default});
  }
}

1;

