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
    'pipeline_name'     => 'DemultiplexIlluminaFastq',
    'seqrun_source'     => undef,
    'seqrun_local_dir'  => undef,
    'seqrun_server'     => undef,
    'base_work_dir'     => undef,
    'seqrun_user'       => undef,
    'checksum_type'     => 'md5',
    'read_offset'       => 1,
    'index_offset'      => 0,
    'bcl2fastq_exe'     => undef,
    'bcl2fastq_options' => '{"-r":1,"-w":1,"-p":1,"--barcode-mismatches":1,"--auto-set-to-zero-barcode-mismatches":"" }',
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
      'seqrun_server' => $self->o('seqrun_server'),
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
    -analysis_capacity => 10,
    -parameters => {
      'seqrun_source'    => $self->o('seqrun_source'),
      'seqrun_user'      => $self->o('seqrun_user'),
      'seqrun_server'    => $self->o('seqrun_server'),
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
    -flow_into => {
      1 => ['find_project_factory'],
    },
  };
  
  push @pipeline, {
    -logic_name  => 'find_project_factory',
    -module      => 'ehive.runnable.jobfactory.SampleSheetProjectFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -flow_into   => {
          '2->A' => ['find_sample_index_length_factory'],
          'A->1' => ['mark_seqrun_status_in_seed_table'],
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
          '2->A' => ['calculate_bases_mask'],
          'A->1' => ['validate_all_lanes_for_project'],
      },
  };

  push @pipeline, {
      -logic_name   => 'calculate_bases_mask',
      -module       => 'ehive.runnable.process.CalculateBasesMask',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -parameters  => {
        'seqrun_local_dir' => $self->o('seqrun_local_dir'),
        'read_offset'      => $self->o('read_offset'),
        'index_offset'     => $self->o('index_offset'),
        },
      -flow_into    => {
          1 => ['run_bcl2fastq']
      },
  };

  push @pipeline, {
      -logic_name   => 'run_bcl2fastq',
      -module       => 'ehive.runnable.process.RunBcl2Fastq',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -parameters  => {
        'seqrun_local_dir'  => $self->o('seqrun_local_dir'),
        'base_work_dir'     => $self->o('base_work_dir'),
        'base_fastq_dir'    => $self->o('base_fastq_dir'),
        'bcl2fastq_exe'     => $self->o('bcl2fastq_exe'),
        'bcl2fastq_options' => $self->o('bcl2fastq_options'),
        },
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
        'seqrun_local_dir' => $self->o('seqrun_local_dir'),
        },
      -flow_into => {
          1 => [ '?accu_name=project_fastq&accu_address={fastq_dir}&accu_input_variable=barcode_qc_stats' ],
        },
  };

  push @pipeline, {
      -logic_name   => 'validate_all_lanes_for_project',
      -module       => 'ehive.runnable.process.ValidateAllLanesForProject',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
  };
  
  push @pipeline, {
      -logic_name   => 'mark_seqrun_status_in_seed_table',
      -module       => 'ehive.runnable.process.ChangePipelineSeedStatus',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -parameters  => {
        'status'        => 'COMPLETED',
        'pipeline_name' => $self->o('pipeline_name'),
        },
  };

  
  return \@pipeline;
}

1;