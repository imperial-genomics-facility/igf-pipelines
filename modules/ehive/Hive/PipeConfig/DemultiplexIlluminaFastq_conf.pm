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
      '2->A' => ['seqrun_transfer_factory'],
      'A->1' => ['check_seqrun_file'],
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
    -flow_into => {
      1 => [ '?accu_name=seqrun_file&accu_address={seqrun_igf_id}&accu_input_variable=seqrun_file_name' ],
    },
  };
  push @pipeline, {
    -logic_name  => 'check_seqrun_file',
    -module      => 'ehive.runnable.IGFBaseProcess',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -parameters => {
    }, 
    -flow_into => {
      1 => ['check_samplesheet'],
    },
  };
  push @pipeline, {
    -logic_name  => 'check_samplesheet',
    -module      => 'ehive.runnable.process.CheckAndProcessSampleSheet',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -parameters  => {
        'seqrun_local_dir' => $self->o{'seqrun_local_dir'},
        'base_work_dir' => $self->o('base_work_dir'),
    },
  -flow_into => {
      1 => ['find_project_factory'],
    },
  };
  push @pipeline, {
      -logic_name  => 'find_project_factory',
      -module      => 'ehive.runnable.jobfactory.SampleSheetProjectFactory',
      -language    => 'python3',
      -meadow_type => 'LOCAL',
      -flow_into =>{
          2 => ['find_flowcell_lane_factory']
      },
  };
  push @pipeline, {
      -logic_name  => 'find_flowcell_lane_factory',
      -module      => 'ehive.runnable.jobfactory.SampleSheetFlowcellFactory',
      -language    => 'python3',
      -meadow_type => 'LOCAL',
      -flow_into   => {
          2 => ['find_sample_index_length_factory'],
      },
  };
  push @pipeline, {
      -logic_name  => 'find_sample_index_length_factory',
      -module      => 'ehive.runnable.jobfactory.SamplesheetFilterAndIndexFactory',
      -language    => 'python3',
      -meadow_type => 'LOCAL',
      -parameters  => {
        'base_work_dir' => $self->o('base_work_dir'),
        },
      -flow_into   => {
          2 => ['calculate_bases_mask'],
      },
  };
  push @pipeline, {
      -logic_name   => 'calculate_bases_mask',
      -module       => 'ehive.runnable.process.CalculateBasesMask',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -flow_into    => {
          1 => ['run_bcl2fastq']
      },
  };
  push @pipeline, {
      -logic_name   => 'run_bcl2fastq',
      -module       => 'ehive.runnable.process.RunBcl2Fastq',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -flow_into    => {
          1 => ['check_demultiplexing_barcode'],
      },
  };
  push @pipeline, {
      -logic_name   => 'check_demultiplexing_barcode',
      -module       => 'ehive.runnable.process.CheckIndexStats',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -parameters  => {
        'base_work_dir' => $self->o('base_work_dir'),
        },
      -flow_into => {
      1 => [ '?accu_name=fastq_files&accu_address={fastq_dir}&accu_input_variable=barcode_qc_stats' ],
        },
  };
};
  
 
  return \@pipeline
}

1;
