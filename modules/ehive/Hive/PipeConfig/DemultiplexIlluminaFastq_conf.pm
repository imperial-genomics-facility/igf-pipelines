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
    'pipeline_name'                  => 'DemultiplexIlluminaFastq',
    'demultiplexing_pipeline_name'   => 'DemultiplexIlluminaFastq',
    'seqrun_source'                  => undef,
    'seqrun_local_dir'               => undef,
    'seqrun_server'                  => undef,
    'base_work_dir'                  => undef,
    'base_results_dir'               => undef,
    'seqrun_user'                    => undef,
    'template_dir'                   => undef,
    'checksum_type'                  => 'md5',
    'read_offset'                    => 1,
    'index_offset'                   => 0,
    'bcl2fastq_exe'                  => undef,
    'bcl2fastq_options'              => '{"-r" : "8","-w" : "4","-p" : "8","--barcode-mismatches" : "1","--auto-set-to-zero-barcode-mismatches":"","--create-fastq-for-index-reads":""}',
    'singlecell_options'             => '{"--minimum-trimmed-read-length=8":"","--mask-short-adapter-reads=8":""}',
    'reset_mask_short_adapter_reads' => 1,
    'fastqc_exe'                     => undef,
    'fastqc_options'                 => '{"-q" : "","--noextract" : "","-f" : "fastq","-k" : "7","-t" : "1"}',
    'irods_exe_dir'                  => undef,
    'fastqscreen_exe'                => undef,
    'fastqscreen_options'            => '{"--aligner" : "bowtie2","--force" : "","--quiet" : "","--subset" : "100000","--threads" : "1"}',
    'fastqscreen_conf'               => undef,
    'multiqc_exe'                    => undef,
    'tool_order_list'                => ['bcl2fastq','fastqc','fastq_screen'],
    'multiqc_template'               => undef,
    'user_info_file'                 => undef,
    'multiqc_options'                => '{"--zip-data-dir" : ""}',
    'cleanup_bcl_dir'                => 0,
    'singlecell_tag'                 => '10X',
    'seqruninfofile'                 => 'seqruninfofile.json',
    'samplereadcountfile'            => 'samplereadcountfile.json',
    'custom_bases_mask'              => undef,
    'use_ephemeral_space'            => 0,
  };
}


sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},            # here we inherit anything from the base class
        'singlecell_tag'     => $self->o('singlecell_tag'),
    };
}


sub pipeline_analyses {
  my ($self) = @_;
  my @pipeline;
  
  push @pipeline, {
    -logic_name          => 'find_new_sequencing_runs',
    -module              => 'ehive.runnable.jobfactory.PipeseedFactory',
    -language            => 'python3',
    -meadow_type         => 'LOCAL',
    -parameters          => {
      'pipeline_name' => $self->o('pipeline_name'),
    },
    -flow_into           => {
      2 => ['seqrun_transfer_factory'],
    },
  };

  push @pipeline, {
    -logic_name          => 'seqrun_transfer_factory',
    -module              => 'ehive.runnable.jobfactory.SeqrunFileFactory',
    -language            => 'python3',
    -meadow_type         => 'PBSPro',
    -rc_name             => '1Gb',
    -analysis_capacity   => 2,
    -parameters          => {
      'seqrun_source' => $self->o('seqrun_source'),
      'seqrun_user'   => $self->o('seqrun_user'),
      'seqrun_server' => $self->o('seqrun_server'),
    },
    -flow_into           => {
      '2->A' => ['transfer_seqrun_file'],
      'A->1' => ['check_samplesheet'],
    },
  };

  push @pipeline, {
    -logic_name          => 'transfer_seqrun_file',
    -module              => 'ehive.runnable.process.demultiplexing.TransferAndCheckRemoteBclFile',
    -language            => 'python3',
    -meadow_type         => 'PBSPro',
    -rc_name             => '1Gb',
    -analysis_capacity   => 10,
    -parameters          => {
      'seqrun_source'      => $self->o('seqrun_source'),
      'seqrun_user'        => $self->o('seqrun_user'),
      'seqrun_server'      => $self->o('seqrun_server'),
      'seqrun_local_dir'   => $self->o('seqrun_local_dir'),
      'checksum_type'      => $self->o('checksum_type'),
    },
  };

  push @pipeline, {
    -logic_name          => 'check_samplesheet',
    -module              => 'ehive.runnable.process.demultiplexing.CheckAndProcessSampleSheet',
    -language            => 'python3',
    -meadow_type         => 'PBSPro',
    -rc_name             => '1Gb',
    -analysis_capacity   => 2,
    -parameters          => {
        'seqrun_local_dir'  => $self->o('seqrun_local_dir'),
        'base_work_dir'     => $self->o('base_work_dir'),
        'singlecell_tag'    => $self->o('singlecell_tag'),
    },
    -flow_into           => {
        1 => WHEN('#project_type# eq #singlecell_tag#' => ['change_single_cell_barcodes'],
	                ELSE ['find_project_factory'],),
    },
  };

  push @pipeline, {
    -logic_name          => 'change_single_cell_barcodes',
    -module              => 'ehive.runnable.process.demultiplexing.ReplaceSingleCellBarcodes',
    -language            => 'python3',
    -meadow_type         => 'PBSPro',
    -rc_name             => '1Gb',
    -analysis_capacity   => 2,
    -parameters          => {
         'base_work_dir'            => $self->o('base_work_dir'),
         'singlecell_tag'           => $self->o('singlecell_tag'),
         'single_cell_barcode_file' => $self->o('single_cell_barcode_file'),
    },
    -flow_into           => {
      1 => ['find_project_factory'],
    },
  };

  push @pipeline, {
    -logic_name          => 'find_project_factory',
    -module              => 'ehive.runnable.jobfactory.SampleSheetProjectFactory',
    -language            => 'python3',
    -meadow_type         => 'PBSPro',
    -rc_name             => '1Gb',
    -analysis_capacity   => 2,
    -flow_into           => {
          '2->A' => ['find_sample_index_length_factory'],
          'A->1' => ['mark_seqrun_status_in_seed_table'],
      },
  };

  push @pipeline, {
      -logic_name        => 'find_sample_index_length_factory',
      -module            => 'ehive.runnable.jobfactory.SamplesheetFilterAndIndexFactory',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '1Gb',
      -analysis_capacity => 2,
      -parameters        => {
        'base_work_dir' => $self->o('base_work_dir'),
        },
      -flow_into         => {
          '2->A' => ['calculate_bases_mask'],
          'A->1' => ['validate_all_lanes_for_project'],
      },
  };

  push @pipeline, {
      -logic_name        => 'calculate_bases_mask',
      -module            => 'ehive.runnable.process.demultiplexing.CalculateBasesMask',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '1Gb',
      -analysis_capacity => 2,
      -parameters        => {
        'seqrun_local_dir' => $self->o('seqrun_local_dir'),
        'read_offset'      => $self->o('read_offset'),
        'index_offset'     => $self->o('index_offset'),
        'custom_bases_mask' => $self->o('custom_bases_mask'),
        },
      -flow_into         => {
          1 => ['run_bcl2fastq']
      },
  };

  push @pipeline, {
      -logic_name        => 'run_bcl2fastq',
      -module            => 'ehive.runnable.process.demultiplexing.RunBcl2Fastq',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '32Gb16t',
      -analysis_capacity => 8,
      -parameters        => {
        'seqrun_local_dir'  => $self->o('seqrun_local_dir'),
        'base_work_dir'     => $self->o('base_work_dir'),
        'base_fastq_dir'    => $self->o('base_fastq_dir'),
        'bcl2fastq_exe'     => $self->o('bcl2fastq_exe'),
        'bcl2fastq_options' => $self->o('bcl2fastq_options'),
        'singlecell_options'=> $self->o('singlecell_options'),
        'singlecell_tag'    => $self->o('singlecell_tag'),
        'use_ephemeral_space'            => $self->o('use_ephemeral_space'),
        'reset_mask_short_adapter_reads' => $self->o('reset_mask_short_adapter_reads'),
        },
      -flow_into         => {
         1 => WHEN('#bcl2fq_project_type# eq #singlecell_tag#' => ['merge_single_cell_fastq'],
                   ELSE ['check_demultiplexing_barcode'],),
      },
  };

  push @pipeline, {
      -logic_name        => 'merge_single_cell_fastq',
      -module            => 'ehive.runnable.process.demultiplexing.MergeSingleCellFastqFragments',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '4Gb',
      -analysis_capacity => 2,
      -parameters        => {
         'singlecell_tag'    => $self->o('singlecell_tag'),
         'use_ephemeral_space' => $self->o('use_ephemeral_space'),
      },
      -flow_into         => {
           1 => ['check_demultiplexing_barcode']
      },
  };

  push @pipeline, {
      -logic_name        => 'check_demultiplexing_barcode',
      -module            => 'ehive.runnable.process.demultiplexing.CheckIndexStats',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '4Gb',
      -analysis_capacity => 8,
      -parameters        => {
        'seqrun_local_dir' => $self->o('seqrun_local_dir'),
        'use_ephemeral_space' => $self->o('use_ephemeral_space'),
        },
      -flow_into         => {
          1 => [ '?accu_name=project_fastq&accu_address={fastq_dir}&accu_input_variable=barcode_qc_stats' ],
        },
  };

  push @pipeline, {
      -logic_name        => 'validate_all_lanes_for_project',
      -module            => 'ehive.runnable.process.demultiplexing.ValidateAllLanesForProject',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '4Gb',
      -analysis_capacity => 2,
      -parameters        => {
        'project_fastq'     => '#project_fastq#',
        },
      -flow_into         => {
          1 => WHEN(
                 '#project_status# eq "PASS"' => [ 'create_remote_project_access' ],
           ),
        },
  };

  push @pipeline, {
      -logic_name        => 'create_remote_project_access',
      -module            => 'ehive.runnable.process.CreateRemoteAccessForProject',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '4Gb',
      -parameters        => {
        'template_dir'        => $self->o('template_dir'),
        'remote_project_path' => $self->o('remote_project_path'),
        'remote_host'         => $self->o('remote_host'),
        'remote_user'         => $self->o('seqrun_user'),
        'seqruninfofile'      => $self->o('seqruninfofile'),
        'samplereadcountfile' => $self->o('samplereadcountfile'),
        'use_ephemeral_space' => $self->o('use_ephemeral_space'),
        },
      -flow_into         => {
          1 => ['project_fastqdir_factory'],
      },
  };

  push @pipeline, {
      -logic_name        => 'project_fastqdir_factory',
      -module            => 'ehive.runnable.jobfactory.ProjectFastqdirFactory',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '4Gb',
      -analysis_capacity => 2,
      -flow_into         => {
          '2->A' => ['collect_fastq_to_db_collection' ],
          'A->1' => ['prepare_and_copy_qc_page_for_project'],
          '2->B' => ['undetermined_fastq_factory'],
          'B->1' => ['prepare_and_copy_qc_page_for_undetermined'],
      },
  };

  push @pipeline, {
      -logic_name        => 'collect_fastq_to_db_collection',
      -module            => 'ehive.runnable.process.demultiplexing.CollectFastqToDbCollection',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '32Gb4t',
      -analysis_capacity => 8,
      -parameters        => {
         'singlecell_tag'    => $self->o('singlecell_tag'),
       },
      -flow_into         => {
          1 => ['upload_fastq_dir_to_irods']
      },
  };

  push @pipeline, {
      -logic_name        => 'upload_fastq_dir_to_irods',
      -module            => 'ehive.runnable.process.demultiplexing.UploadFastqToIrods',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 4,
      -parameters        => {
        'irods_exe_dir' => $self->o('irods_exe_dir'),
        'use_ephemeral_space' => $self->o('use_ephemeral_space'),
        },
      -flow_into         => {
          1 => ['known_fastq_factory']
      },
  };

  push @pipeline, {
      -logic_name        => 'known_fastq_factory',
      -module            => 'ehive.runnable.jobfactory.FastqFileFactory',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 2,
      -parameters        => {
        'filter_keyword'   => 'Undetermined*',
        'required_keyword' => undef,
        },
      -flow_into         => {
          '2->A' => ['run_fastqc_for_known_fastq'],
          'A->1' => ['collect_qc_data_for_known_fastq'],
      },
  };

  push @pipeline, {
      -logic_name        => 'run_fastqc_for_known_fastq',
      -module            => 'ehive.runnable.process.demultiplexing.RunFastqc',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 30,
      -parameters        => {
        'base_results_dir' => $self->o('base_results_dir'),
        'fastqc_exe'       => $self->o('fastqc_exe'),
        'fastqc_options'   => $self->o('fastqc_options'),
        'tag'              => 'known',
        'use_ephemeral_space' => $self->o('use_ephemeral_space'),
        },
      -flow_into         => {
          1 => {'copy_fastqc_results_to_remote' =>
                  {'file'            => '#fastqc_html#',
                   'lane_index_info' => '#lane_index_info#',
                   'sample_name'     => '#sample_name#',
                   'fastqc'          => '#fastqc#'},
               },
      },
  };

  push @pipeline, {
      -logic_name        => 'copy_fastqc_results_to_remote',
      -module            => 'ehive.runnable.process.demultiplexing.CopyQCFileToRemote',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 2,
      -parameters        => {
        'file'                => '#file#',
        'tag'                 => 'known',
        'analysis_label'      => 'fastqc',
        'dir_label'           => '#lane_index_info#',
        'sample_label'        => '#sample_name#',
        'remote_user'         => $self->o('seqrun_user'),
        'remote_host'         => $self->o('remote_host'),
        'remote_project_path' => $self->o('remote_project_path'),
        },
      -flow_into         => {
          1 => {'run_fastqscreen_for_known_fastq' =>
                  {'fastqc_remote_file'=>'#remote_file#'}
               },
      },
  };

  push @pipeline, {
      -logic_name        => 'run_fastqscreen_for_known_fastq',
      -module            => 'ehive.runnable.process.demultiplexing.RunFastqscreen',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 30,
      -parameters        => {
        'base_results_dir'    => $self->o('base_results_dir'),
        'fastqscreen_exe'     => $self->o('fastqscreen_exe'),
        'fastqscreen_options' => $self->o('fastqscreen_options'),
        'fastqscreen_conf'    => $self->o('fastqscreen_conf'),
        'tag'                 => 'known',
        'use_ephemeral_space' => $self->o('use_ephemeral_space'),
        },
      -flow_into         => {
          1 => {'copy_fastqscreen_results_to_remote' => 
                   {'file'        => '#fastqscreen_html#',
                    'fastqscreen' => '#fastqscreen#'}},
      },
  };

  push @pipeline, {
      -logic_name        => 'copy_fastqscreen_results_to_remote',
      -module            => 'ehive.runnable.process.demultiplexing.CopyQCFileToRemote',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 2,
      -parameters        => {
        'file'                => '#file#',
        'tag'                 => 'known',
        'analysis_label'      => 'fastqscreen',
        'dir_label'           => '#lane_index_info#',
        'sample_label'        => '#sample_name#',
        'remote_user'         => $self->o('seqrun_user'),
        'remote_host'         => $self->o('remote_host'),
        'remote_project_path' => $self->o('remote_project_path'),
        },
      -flow_into         => {
          1 => ['?accu_name=known_fastqc&accu_address={fastq_file}&accu_input_variable=fastqc',
                '?accu_name=known_fastscreen&accu_address={fastq_file}&accu_input_variable=fastqscreen',
                '?accu_name=known_remote_fastqc&accu_address={fastq_file}&accu_input_variable=fastqc_remote_file',
                '?accu_name=known_remote_fastscreen&accu_address={fastq_file}&accu_input_variable=remote_file',
               ],
      },
  };

  push @pipeline, {
      -logic_name        => 'collect_qc_data_for_known_fastq',
      -module            => 'ehive.runnable.process.demultiplexing.CollectQcForFastqDir',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 2,
      -parameters        => {
          'fastqc_info'        => '#known_fastqc#',
          'fastqscreen_info'   => '#known_fastscreen#',
          'remote_fastqc_info' => '#known_remote_fastqc#',
          'remote_fastqs_info' => '#known_remote_fastscreen#',
      },
      -flow_into         => {
          1 => ['run_multiqc_for_know_fastq'],
      },
  };

  push @pipeline, {
      -logic_name        => 'run_multiqc_for_know_fastq',
      -module            => 'ehive.runnable.process.demultiplexing.RunMutiQC',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 8,
      -parameters        => {
        'qc_files'         => '#qc_outputs#',
        'base_results_dir' => $self->o('base_results_dir'),
        'multiqc_exe'      => $self->o('multiqc_exe'),
        'multiqc_options'  => $self->o('multiqc_options'),
        'tag'              => 'known',
        'tool_order_list'  => $self->o('tool_order_list'),
        'multiqc_template_file' => $self->o('multiqc_template'),
        'use_ephemeral_space'   => $self->o('use_ephemeral_space'),
        },
      -flow_into         => {
          1 => {'copy_known_multiqc_to_remote' =>
                  {'file'            => '#multiqc_html#',
                   'lane_index_info' => '#lane_index_info#',
                   'qc_outputs'      => '#qc_outputs#'},
               },
      },
  };

  push @pipeline, {
      -logic_name        => 'copy_known_multiqc_to_remote',
      -module            => 'ehive.runnable.process.demultiplexing.CopyQCFileToRemote',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 2,
      -parameters        => {
        'tag'                 => 'known',
        'file'                => '#file#',
        'analysis_label'      => 'multiqc',
        'dir_label'           => '#lane_index_info#',
        'remote_user'         => $self->o('seqrun_user'),
        'remote_host'         => $self->o('remote_host'),
        'remote_project_path' => $self->o('remote_project_path'),
        },
      -flow_into         => {
          1 => {'prepare_and_copy_qc_page_for_lane' =>
                   {'multiqc_remote_file' => '#remote_file#',
                    'qc_files'            => '#qc_outputs#', },
               },
      },
  };

  push @pipeline, {
      -logic_name        => 'prepare_and_copy_qc_page_for_lane',
      -module            => 'ehive.runnable.process.demultiplexing.PrepareQcPageForRemote',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 4,
      -parameters        => {
        'qc_files'            => '#qc_files#',
        'template_dir'        => $self->o('template_dir'),
        'remote_host'         => $self->o('remote_host'),
        'remote_user'         => $self->o('seqrun_user'),
        'page_type'           => 'sample',
        'remote_project_path' => $self->o('remote_project_path'),
        'singlecell_tag'      => $self->o('singlecell_tag'),
        'use_ephemeral_space' => $self->o('use_ephemeral_space'),
        },
      -flow_into         => {
          1 => ['?accu_name=laneqc_known&accu_address={fastq_dir}&accu_input_variable=qc_file_info',],
      },
  };

  push @pipeline, {
      -logic_name        => 'prepare_and_copy_qc_page_for_project',
      -module            => 'ehive.runnable.process.demultiplexing.PrepareQcPageForRemote',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 2,
      -parameters        => {
        'qc_files'            => '#laneqc_known#',
        'template_dir'        => $self->o('template_dir'),
        'remote_host'         => $self->o('remote_host'),
        'remote_user'         => $self->o('seqrun_user'),
        'page_type'           => 'project',
        'remote_project_path' => $self->o('remote_project_path'),
        'use_ephemeral_space' => $self->o('use_ephemeral_space'),
        },
      -flow_into         => {
          1 => ['update_project_info_page'],
      },
  };

  push @pipeline, {
      -logic_name        => 'update_project_info_page',
      -module            => 'ehive.runnable.process.UpdateProjectInfo',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 2,
      -parameters        => {
        'remote_project_path' => $self->o('remote_project_path'),
        'remote_host'         => $self->o('remote_host'),
        'remote_user'         => $self->o('seqrun_user'),
        'seqruninfofile'      => $self->o('seqruninfofile'),
        'samplereadcountfile' => $self->o('samplereadcountfile'),
        'pipeline_name'       => $self->o('demultiplexing_pipeline_name'),
        'use_ephemeral_space' => $self->o('use_ephemeral_space'),
        'analysis_pipeline_name' => $self->o('analysis_pipeline_name'),
      },
      -flow_into         => {
          1 => ['send_email_notification'],
      },
  };

  push @pipeline, {
      -logic_name        => 'send_email_notification',
      -module            => 'ehive.runnable.process.demultiplexing.SendEmailToUser',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 2,
      -parameters        => {
        'template_dir'    => $self->o('template_dir'),
        'remote_host'     => $self->o('seqrun_server'),
        'remote_user'     => $self->o('seqrun_user'),
        'user_info_file'  => $self->o('user_info_file'),
        'use_ephemeral_space' => $self->o('use_ephemeral_space'),
        },
  };

  push @pipeline, {
      -logic_name        => 'undetermined_fastq_factory',
      -module            => 'ehive.runnable.jobfactory.FastqFileFactory',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 2,
      -parameters        => {
        'required_keyword' => 'Undetermined*',
        'filter_keyword'   => undef,
        },
      -flow_into         => {
          '2->A' => ['run_fastqc_for_undetermined_fastq'],
          'A->1' => ['collect_qc_data_for_undetermined_fastq'],
      },
  };

  push @pipeline, {
      -logic_name        => 'run_fastqc_for_undetermined_fastq',
      -module            => 'ehive.runnable.process.demultiplexing.RunFastqc',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 10,
      -parameters        => {
        'base_results_dir' => $self->o('base_results_dir'),
        'fastqc_exe'       => $self->o('fastqc_exe'),
        'fastqc_options'   => $self->o('fastqc_options'),
        'tag'              => 'undetermined',
        'use_ephemeral_space' => $self->o('use_ephemeral_space'),
        },
      -flow_into         => {
          1 => ['run_fastqscreen_for_undetermined_fastq']
      },
  };

  push @pipeline, {
      -logic_name        => 'run_fastqscreen_for_undetermined_fastq',
      -module            => 'ehive.runnable.process.demultiplexing.RunFastqscreen',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 10,
      -parameters        => {
        'base_results_dir'    => $self->o('base_results_dir'),
        'fastqscreen_exe'     => $self->o('fastqscreen_exe'),
        'fastqscreen_options' => $self->o('fastqscreen_options'),
        'fastqscreen_conf'    => $self->o('fastqscreen_conf'),
        'tag'                 => 'undetermined',
        'use_ephemeral_space' => $self->o('use_ephemeral_space'),
        },
      -flow_into         => {
           1 => ['?accu_name=undetermined_fastqc&accu_address={fastq_file}&accu_input_variable=fastqc',
                 '?accu_name=undetermined_fastscreen&accu_address={fastq_file}&accu_input_variable=fastqscreen',
                ],
      },
  };

  push @pipeline, {
      -logic_name        => 'collect_qc_data_for_undetermined_fastq',
      -module            => 'ehive.runnable.process.demultiplexing.CollectQcForFastqDir',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -parameters        => {
          'fastqc_info'      => '#undetermined_fastqc#',
          'fastqscreen_info' => '#undetermined_fastscreen#',
      },
      -flow_into         => {
          1 => ['run_multiqc_for_undetermined_fastq'],
      },
  };

  push @pipeline, {
      -logic_name        => 'run_multiqc_for_undetermined_fastq',
      -module            => 'ehive.runnable.process.demultiplexing.RunMutiQC',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 8,
      -parameters        => {
        'qc_files'         => '#qc_outputs#',
        'base_results_dir' => $self->o('base_results_dir'),
        'multiqc_exe'      => $self->o('multiqc_exe'),
        'multiqc_options'  => $self->o('multiqc_options'),
        'tag'              => 'undetermined',
        'tool_order_list'  => $self->o('tool_order_list'),
        'multiqc_template_file' => $self->o('multiqc_template'),
        'use_ephemeral_space'   => $self->o('use_ephemeral_space'),
        },
      -flow_into         => {
          1 => {'copy_undetermined_multiqc_to_remote' =>
                  {'file'            => '#multiqc_html#',
                   'lane_index_info' => '#lane_index_info#'},
               },
      },
  };

  push @pipeline, {
      -logic_name        => 'copy_undetermined_multiqc_to_remote',
      -module            => 'ehive.runnable.process.demultiplexing.CopyQCFileToRemote',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 2,
      -parameters        => {
        'file'                => '#file#',
        'tag'                 => 'undetermined',
        'analysis_label'      => 'multiqc',
        'dir_label'           => '#lane_index_info#',
        'remote_user'         => $self->o('seqrun_user'),
        'remote_host'         => $self->o('remote_host'),
        'remote_project_path' => $self->o('remote_project_path'),
        },
      -flow_into         => {
          1 => ['?accu_name=multiqc_undetermined&accu_address={fastq_dir}&accu_input_variable=remote_file',],
      },
  };

  push @pipeline, {
      -logic_name        => 'prepare_and_copy_qc_page_for_undetermined',
      -module            => 'ehive.runnable.process.demultiplexing.PrepareQcPageForRemote',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -analysis_capacity => 2,
      -parameters        => {
        'qc_files'            => '#multiqc_undetermined#',
        'template_dir'        => $self->o('template_dir'),
        'remote_host'         => $self->o('remote_host'),
        'remote_user'         => $self->o('seqrun_user'),
        'page_type'           => 'undetermined',
        'remote_project_path' => $self->o('remote_project_path'),
        'fastq_dir'           => undef,
        'use_ephemeral_space' => $self->o('use_ephemeral_space'),
        },
  };

  push @pipeline, {
      -logic_name        => 'mark_seqrun_status_in_seed_table',
      -module            => 'ehive.runnable.process.ChangePipelineSeedStatus',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -parameters        => {
        'igf_id'        => '#seqrun_igf_id#',
        'task_id'       => '#seqrun_igf_id#',
        'new_status'    => 'FINISHED',
        'pipeline_name' => $self->o('pipeline_name'),
        },
      -flow_into         => {
        1 => ['remove_bcl_directory'],  
      },
  };

  push @pipeline, {
      -logic_name        => 'remove_bcl_directory',
      -module            => 'ehive.runnable.process.CleanupDirOrFile',
      -language          => 'python3',
      -meadow_type       => 'PBSPro',
      -rc_name           => '8Gb2t',
      -parameters        => {
        'seqrun_local_dir' => $self->o('seqrun_local_dir'),
        'path'             => '#seqrun_local_dir#/#seqrun_igf_id#',
        'cleanup_status'   => $self->o('cleanup_bcl_dir'),
        },
  };

  return \@pipeline;
}

1;
