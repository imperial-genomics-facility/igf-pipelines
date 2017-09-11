=head1 NAME
    ehive::Hive::PipeConfig::DemultiplexIlluminaFastq_conf

=cut

package ehive::Hive::PipeConfig::DemultiplexIlluminaFastq_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('ehive::Hive::PipeConfig::IGFBasePipe_conf');

sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },
    'pipeline_name' => 'DemultiplexIlluminaFastq',
    'seed_id_label' => undef,
    'seqrun_id_label' => undef,
    'seqrun_source' => undef,
    'seqrun_server' => undef,
    'seqrun_md5_type' => undef,
    'seqrun_local_dir' => undef,
    'checksum_type' => undef,
  };
}

sub pipeline_analyses {
  my ($self) = @_;
  my @pipeline;
  
  push @pipeline, {
    -logic_name  => 'find_new_sequencing_runs',
    -module      => 'ehive.runnable.jobfactory.PipeseedFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -parameters  =>{
      'pipeline_name' => $self->o('pipeline_name'),
      'seed_id_label' => $self->o('seed_id_label'),
      'seqrun_id_label'=> $self->o('seqrun_id_label'),
    },
    -flow_into => {
      2 => ['seqrun_transfer_factory'],    
    },
  };
  push @pipeline, {
    -logic_name  => 'seqrun_transfer_factory',
    -module      => 'ehive.runnable.jobfactory.SeqrunFileFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -parameters  => {
      'seqrun_source' => $self->o('seqrun_source'),
      'seqrun_server' => $self->o('seqrun_server'),
      'seqrun_user'   => $self->o('seqrun_user'),
      'seqrun_md5_type' => $self->o('seqrun_md5_type'),
    },
    -flow_into => {
      2 => ['transfer_seqrun_file'],    
    },
  };
  push @pipeline, {
    -logic_name  => 'transfer_seqrun_file',
    -module      => 'ehive.runnable.proccess.TransferAndCheckRemoteFile',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -parameters => {
      'seqrun_source' => $self->o('seqrun_source'),
      'seqrun_server' => $self->o('seqrun_server'),
      'seqrun_local_dir' => $self->o('seqrun_local_dir'),
      'checksum_type' => $self->o('checksum_type'),
    },
  };
  
 
  return \@pipeline
}

1;
