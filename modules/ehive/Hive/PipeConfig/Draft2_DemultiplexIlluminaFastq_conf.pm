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
    'base_work_dir' => undef,
    'seqrun_user' => undef,
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
    -module      => 'ehive.runnable.proccess.TransferAndCheckRemoteBclFile',
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
          '2->A' => ['find_flowcell_lane_factory'],
          'A->1' => ['mark_seqrun_status_in_seed_table'],
      },
  };
  push @pipeline, {
      -logic_name  => 'find_flowcell_lane_factory',
      -module      => 'ehive.runnable.jobfactory.SampleSheetFlowcellFactory',
      -language    => 'python3',
      -meadow_type => 'LOCAL',
      -flow_into   => {
          '2->A' => ['find_sample_index_length_factory'],
          'A->1' => ['validate_all_lanes_for_project'],
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
          1 => [ '?accu_name=project_fastq&accu_address={fastq_dir}&accu_input_variable=barcode_qc_stats' ],
        },
  };
  push @pipeline, {
      -logic_name   => 'validate_all_lanes_for_project',
      -module       => 'ehive.runnable.process.ValidateAllLanesForProject',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -flow_into   => {
          1 => WHEN( 
                 '#project_fastq#' => { 'project_fastqdir_factory' => INPUT_PLUS(),},
                 ELSE { 'mark_seqrun_failed' => INPUT_PLUS(),},
           ),
      },
  };

  push @pipeline, {
      -logic_name   => 'mark_seqrun_failed',
      -module       => 'ehive.runnable.process.ChangePipelineSeedStatus',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -parameters  => {
        'status' => 'FAILED*',
        },
  };

  push @pipeline, {
      -logic_name   => 'project_fastqdir_factory',
      -module       => 'ehive.runnable.process.ProjectFastqdirFactory',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -flow_into    => {
          '2->A' => ['collect_fastq_to_db_collection' ],
          'A->1' => ['run_multiqc_for_know_fastq'],
          '2->B' => ['undetermined_fastq_factory'],
          'B->1' => ['run_multiqc_for_undetermined_fastq'],
      },
  };
  
  push @pipeline, {
      -logic_name   => 'collect_fastq_to_db_collection',
      -module       => 'ehive.runnable.process.CollectFastqToDbCollection',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -flow_into    => {
          1 => ['upload_fastq_dir_to_irods']
      },
  };

  push @pipeline, {
      -logic_name   => 'upload_fastq_dir_to_irods',
      -module       => 'ehive.runnable.process.UploadFastqToIrods',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -flow_into    => {
          1 => ['known_fastq_factory']
      },
  };

  push @pipeline, {
      -logic_name   => 'known_fastq_factory',
      -module       => 'ehive.runnable.jobfactory.FastqFileFactory',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -parameters  => {
        'filter_keyword' => 'Undetermined*',
        },
      -flow_into    => {
          '2->A' => ['run_fastqc_for_known_fastq'],
          'A->1' => ['collect_qc_data_for_known_fastq'],
      },
  };

  push @pipeline, {
      -logic_name   => 'run_fastqc_for_known_fastq',
      -module       => 'ehive.runnable.process.RunFastqc',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -flow_into    => {
          1 => ['copy_fastqc_results_to_remote']
      },
  };

  push @pipeline, {
      -logic_name   => 'copy_fastqc_results_to_remote',
      -module       => 'ehive.runnable.process.CopyQCFileToRemote',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -parameters  => {
        'tag' => 'known',
        'label'=>'fastqc',
        },
      -flow_into    => {
          1 => ['run_fastqscreen_for_known_fastq']
      },
  };

  push @pipeline, {
      -logic_name   => 'run_fastqscreen_for_known_fastq',
      -module       => 'ehive.runnable.process.RunFastqscreen',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -flow_into    => {
          1 => ['copy_fastqscreen_results_to_remote']
      },
  };

  push @pipeline, {
      -logic_name   => 'copy_fastqscreen_results_to_remote',
      -module       => 'ehive.runnable.process.CopyQCFileToRemote',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -parameters  => {
        'tag' => 'known',
        'label'=>'fastqscreen',
        },
      -flow_into    => {
          1 => ['?accu_name=known_fastqc&accu_address={fastq_file}&accu_input_variable=fastqc_output',
                '?accu_name=known_fastscreen&accu_address={fastq_file}&accu_input_variable=fastqscreen_output',
               ],
      },
  };

  push @pipeline, {
      -logic_name   => 'collect_qc_data_for_known_fastq',
      -module       => 'ehive.runnable.process.CollectQcForFastqDir',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -flow_into    => {
          1 => ['?accu_name=qc_known&accu_address={fastq_dir}&accu_input_variable=qc_outputs',]
      },
  };

  push @pipeline, {
      -logic_name   => 'run_multiqc_for_know_fastq',
      -module       => 'ehive.runnable.process.RunMutiQC',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -flow_into    => {
          1 => ['copy_known_multiqc_to_remote']
      },
  };

  push @pipeline, {
      -logic_name   => 'copy_known_multiqc_to_remote',
      -module       => 'ehive.runnable.process.CopyQCFileToRemote',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -parameters  => {
        'tag' => 'known',
        'label'=>'multiqc',
        },
  };

  push @pipeline, {
      -logic_name   => 'undetermined_fastq_factory',
      -module       => 'ehive.runnable.jobfactory.FastqFileFactory',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -parameters  => {
        'required_keyword' => 'Undetermined*',
        },
      -flow_into    => {
          '2->A' => ['run_fastqc_for_undetermined_fastq'],
          'A->1' => ['collect_qc_data_for_undetermined_fastq'],
      },
  };

  push @pipeline, {
      -logic_name   => 'run_fastqc_for_undetermined_fastq',
      -module       => 'ehive.runnable.process.RunFastqc',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -flow_into    => {
          1 => ['run_fastqscreen_for_undetermined_fastq']
      },
  };

  push @pipeline, {
      -logic_name   => 'run_fastqscreen_for_undetermined_fastq',
      -module       => 'ehive.runnable.process.RunFastqscreen',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -flow_into    => {
           1 => ['?accu_name=undetermined_fastqc&accu_address={fastq_file}&accu_input_variable=fastqc_output',
                 '?accu_name=undetermined_fastscreen&accu_address={fastq_file}&accu_input_variable=fastqscreen_output',
                ],
      },
  };

  push @pipeline, {
      -logic_name   => 'collect_qc_data_for_undetermined_fastq',
      -module       => 'ehive.runnable.process.CollectQcForFastqDir',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -flow_into    => {
          1 => ['?accu_name=qc_undetermined&accu_address={fastq_dir}&accu_input_variable=qc_outputs',]
      },
  };

  push @pipeline, {
      -logic_name   => 'run_multiqc_for_undetermined_fastq',
      -module       => 'ehive.runnable.process.RunMutiQC',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -flow_into    => {
          1 => ['copy_undetermined_multiqc_to_remote']
      },
  };

  push @pipeline, {
      -logic_name   => 'copy_undetermined_multiqc_to_remote',
      -module       => 'ehive.runnable.process.CopyQCFileToRemote',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -parameters  => {
        'tag' => 'undetermined',
        'label'=>'multiqc',
        },
  };

  push @pipeline, {
      -logic_name   => 'mark_seqrun_status_in_seed_table',
      -module       => 'ehive.runnable.process.ChangePipelineSeedStatus',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -parameters  => {
        'status' => 'COMPLETED',
        },
  };

  return \@pipeline;
}

1;