=head1 NAME
    ehive::Hive::PipeConfig::TestDemultiplexIlluminaFastq_conf
=cut

package ehive::Hive::PipeConfig::TestDemultiplexIlluminaFastq_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('ehive::Hive::PipeConfig::IGFBasePipe_conf');

sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },
    'pipeline_name'  => 'DemultiplexIlluminaFastq',
    'seqrun_source'  => undef,
    'seqrun_local_dir' => undef,
    'base_work_dir' => undef,
    'seqrun_user'   => undef,
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
      'seqrun_user'   => $self->o('seqrun_user'),
    },
    -flow_into => {
      '2->A' => ['transfer_seqrun_file'],
      'A->1' => ['check_samplesheet'],
    },
  };

  push @pipeline, {
    -logic_name  => 'transfer_seqrun_file',
    -module      => 'ehive.runnable.process.TransferAndCheckRemoteBclFile',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '500Mb',
    -parameters => {
      'seqrun_source'    => $self->o('seqrun_source'),
      'seqrun_user'      => $self->o('seqrun_user'),
      'seqrun_local_dir' => $self->o('seqrun_local_dir'),
      'checksum_type'    => $self->o('checksum_type'),
    },
    -flow_into => {
      1 => [ '?accu_name=seqrun_file&accu_address={seqrun_igf_id}&accu_input_variable=seqrun_file_name' ],
    },
  };

  push @pipeline, {
    -logic_name  => 'check_samplesheet',
    -module      => 'ehive.runnable.process.CheckAndProcessSampleSheet',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -parameters  => {
        'seqrun_local_dir' => $self->o('seqrun_local_dir'),
        'base_work_dir' => $self->o('base_work_dir'),
    },
  };
  return \@pipeline;
}

1;