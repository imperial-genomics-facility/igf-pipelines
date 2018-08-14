=head1 NAME
    ehive::Hive::PipeConfig::PrimaryAnalysis_conf
=cut


package ehive::Hive::PipeConfig::PrimaryAnalysis_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('ehive::Hive::PipeConfig::IGFBasePipe_conf');

sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },                                       # here we inherit anything from the base class
    'pipeline_name'       => 'PrimaryAnalysis',
    'pipeseed_mode'       => 'alignment',
    'singlecell_source'   => 'TRANSCRIPTOMIC_SINGLE_CELL',
    'tenx_exp_type'       => 'TENX-TRANSCRIPTOME',
    'base_work_dir'       => undef,
    'base_results_dir'    => undef,
    'seqrun_user'         => undef,
    'template_dir'        => undef,
    'checksum_type'       => 'md5',
    'irods_exe_dir'       => undef,
    'cellranger_exe'      => undef,
    'cellranger_param'    => '{"--nopreflight":"","--disable-ui":"","--jobmode":"pbspro","--localcores":"1","--localmem":"1","--mempercore":"4","--maxjobs":"20"}',
    'multiqc_options'     => '{"--zip-data-dir" : ""}',
    'cleanup_bam_dir'     => 0,
    'cram_type'           => 'ANALYSIS_CRAM',
    'multiqc_type'        => 'MULTIQC_HTML',
    'scanpy_type'         => 'SCANPY_RESULTS',
    'samtools_threads'    => 4,
    'cellranger_timeout'  => 43200,
    'java_exe'            => undef,
    'picard_jar'          => undef,
    'java_param'          => '-Xmx4g',
    'copy_input_to_temp'  => 0,
    'multiqc_exe'         => undef,
    'multiqc_options'     => '{"--zip-data-dir" : ""}',
    'reference_fasta_type'        => 'GENOME_FASTA',
    'reference_refFlat'           => 'GENE_REFFLAT',
    'cellranger_collection_table' => 'experiment',
    'cellranger_analysis_name'    => 'cellranger_count',
    'multiqc_analysis_name'       => 'multiqc',
    'scanpy_report_template'      => undef,
  };
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},                              # here we inherit anything from the base class
        'singlecell_source' => $self->o('singlecell_source'),
        'tenx_exp_type' => $self->o('tenx_exp_type'),
    };
}

sub pipeline_analyses {
  my ($self) = @_;
  my @pipeline;
  
  ## collect all experiment seeds
  push @pipeline, {
    -logic_name  => 'find_new_experiment_for_analysis',
    -module      => 'ehive.runnable.jobfactory.PipeseedFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -parameters  => {
      'pipeline_name' => $self->o('pipeline_name'),
      'pipeseed_mode' => $self->o('pipeseed_mode'),
    },
    -flow_into => {
        2 => WHEN('#experiment_type# eq #tenx_exp_type# && #library_source# eq #singlecell_source#' => ['run_cellranger_count_for_experiment'],
                 ELSE ['mark_experiment_finished'],),
    },
  };
  
  ## run cellranger for each experiments
  push @pipeline, {
    -logic_name  => 'run_cellranger_count_for_experiment',
    -module      => 'ehive.runnable.process.alignment.RunCellrangerCount',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 1,
    -parameters  => {
      'cellranger_exe'     => $self->o('cellranger_exe'),
      'cellranger_options' => $self->o('cellranger_param'),
      'base_work_dir'      => $self->o('base_work_dir'),
      'base_results_dir'   => $self->o('base_results_dir'),
      'job_timeout'        => $self->o('cellranger_timeout'),
      },
    -flow_into   => {
        1 => ['load_cellranger_count_output_for_experiment'],
      },
  };
  
  ## load cellranger output for each experiments
  push @pipeline, {
    -logic_name  => 'load_cellranger_count_output_for_experiment',
    -module      => 'ehive.runnable.process.alignment.ProcessCellrangerCountOutput',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 1,
    -parameters  => {
      'base_work_dir'      => $self->o('base_work_dir'),
      'base_results_dir'   => $self->o('base_results_dir'),
      },
    -flow_into   => {
        1 => ['single_cell_analysis_factory_1'],
        #1 => ['upload_cellranger_results_to_irods'],
      },
  };
  
  ## single cell analysis factory 1
  push @pipeline, {
    -logic_name  => 'single_cell_analysis_factory_1',
    -module      => 'ehive.runnable.jobfactory.alignment.AnalysisFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -parameters  => {
      'file_list' => ['#bam_file#'],
      },
    -flow_into   => {
        '2->A' => ['upload_cellranger_results_to_irods',
                   'convert_bam_to_cram',
                   'scanpy_report_generation'
                   ],
        'A->1' => ['picard_aln_summary_for_cellranger'], 
      },
  };
  
  ## upload cellranger resilts to irods
  push @pipeline, {
    -logic_name  => 'upload_cellranger_results_to_irods',
    -module      => 'ehive.runnable.process.alignment.UploadAnalysisResultsToIrods',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'file_list'     => '#analysis_output_list#',
      'irods_exe_dir' => $self->o('irods_exe_dir'),
      'analysis_name' => $self->o('cellranger_analysis_name'),
      'dir_path_list' => ['#sample_igf_id#','#experiment_igf_id#','#analysis_name#'],
      'file_tag'      => '#sample_igf_id#'.' - '.'#experiment_igf_id#'.' - '.'#analysis_name#'.' - '.'#species_name#',
     },
  };
  
  ## convert bam file to cram
  push @pipeline, {
    -logic_name  => 'convert_bam_to_cram',
    -module      => 'ehive.runnable.process.alignment.ConvertBamToCram',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
        'bam_file'        => '#bam_file#',
        'base_result_dir' => $self->o('base_results_dir'),
        'threads'         => $self->o('samtools_threads'),
        'collection_name' => '#experiment_igf_id#',
        'collection_type' => $self->o('cram_type'),
        'collection_table'=> $self->o('cellranger_collection_table'),
        'analysis_name'   => $self->o('cellranger_analysis_name'),
        'tag_name'        => '#species_name#',
        'reference_type'  => $self->o('reference_fasta_type'),
        'copy_input'      => $self->o('copy_input_to_temp'),
     },
     -flow_into   => {
        1 => ['upload_cellranger_cram_to_irods'],
      },
  };
  
  ## upload cram file to irods server
  push @pipeline, {
    -logic_name  => 'upload_cellranger_cram_to_irods',
    -module      => 'ehive.runnable.process.alignment.UploadAnalysisResultsToIrods',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'file_list'     => '#output_cram_list#',
      'irods_exe_dir' => $self->o('irods_exe_dir'),
      'analysis_name' => $self->o('cellranger_analysis_name'),
      'dir_path_list' => ['#sample_igf_id#','#experiment_igf_id#','#analysis_name#'],
      'file_tag'      => '#sample_igf_id#'.' - '.'#experiment_igf_id#'.' - '.'#analysis_name#'.' - '.'#species_name#',
     },
  };
  
  ## scanpy report
  push @pipeline, {
    -logic_name  => 'scanpy_report_generation',
    -module      => 'ehive.runnable.process.alignment.RunScanpy',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'report_template_file'   => $self->o('scanpy_report_template'),
      'base_result_dir'        => $self->o('base_results_dir'),
      'scanpy_collection_type' => $self->o('scanpy_type'),
     },
    -flow_into   => {
        1 => ['copy_scanpy_report_to_remote'],
      },
  };
  
  ## copy scanpy report to remote
  push @pipeline, {
      -logic_name   => 'copy_scanpy_report_to_remote',
      -module       => 'ehive.runnable.process.alignment.CopyAnalysisFilesToRemote',
      -language     => 'python3',
      -meadow_type  => 'PBSPro',
      -rc_name      => '1Gb',
      -analysis_capacity => 2,
      -parameters  => {
        'file_list'           => ['#output_report#'],
        'remote_user'         => $self->o('seqrun_user'),
        'remote_host'         => $self->o('remote_host'),
        'remote_project_path' => $self->o('remote_project_path'),
        },
  };
  
  ## picard alignment summary metrics
  push @pipeline, {
    -logic_name  => 'picard_aln_summary_for_cellranger',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_file'     => '#bam_file#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'CollectAlignmentSummaryMetrics',
      'base_work_dir'  => $self->o('base_work_dir'),
      'copy_input'     => $self->o('copy_input_to_temp'),
      'reference_type' => $self->o('reference_fasta_type'),
     },
    -flow_into   => {
        1 => ['picard_base_dist_summary_for_cellranger'],
      },
  };
  
  ## picard base distribution summary metrics
  push @pipeline, {
    -logic_name  => 'picard_base_dist_summary_for_cellranger',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_file'     => '#bam_file#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'CollectBaseDistributionByCycle',
      'base_work_dir'  => $self->o('base_work_dir'),
      'copy_input'     => $self->o('copy_input_to_temp'),
      'reference_type' => $self->o('reference_fasta_type'),
     },
    -flow_into   => {
        1 => ['picard_gc_bias_summary_for_cellranger'],
      },
  };
  
  ## picard gc bias summary metrics
  push @pipeline, {
    -logic_name  => 'picard_gc_bias_summary_for_cellranger',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_file'     => '#bam_file#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'CollectGcBiasMetrics',
      'base_work_dir'  => $self->o('base_work_dir'),
      'copy_input'     => $self->o('copy_input_to_temp'),
      'reference_type' => $self->o('reference_fasta_type'),
     },
    -flow_into   => {
        1 => ['picard_qual_dist_summary_for_cellranger'],
      },
  };
  
  ## picard quality distribution summary metrics
  push @pipeline, {
    -logic_name  => 'picard_qual_dist_summary_for_cellranger',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_file'     => '#bam_file#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'QualityScoreDistribution',
      'base_work_dir'  => $self->o('base_work_dir'),
      'copy_input'     => $self->o('copy_input_to_temp'),
      'reference_type' => $self->o('reference_fasta_type'),
     },
    -flow_into   => {
        1 => ['picard_rna_metrics_summary_for_cellranger'],
      },
  };
  
  ## picard rna metrics summary metrics
  push @pipeline, {
    -logic_name  => 'picard_rna_metrics_summary_for_cellranger',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_file'     => '#bam_file#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'CollectRnaSeqMetrics',
      'base_work_dir'  => $self->o('base_work_dir'),
      'reference_type' => $self->o('reference_fasta_type'),
      'copy_input'     => $self->o('copy_input_to_temp'),
      'reference_refFlat' => $self->o('reference_refFlat'),
     },
    -flow_into   => {
        1 => ['samtools_flagstat_summary_for_cellranger'],
      },
  };
  
  ## samtools flagstat metrics
  push @pipeline, {
    -logic_name  => 'samtools_flagstat_summary_for_cellranger',
    -module      => 'ehive.runnable.process.alignment.RunSamtools',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'bam_file'         => '#bam_file#',
      'samtools_command' => 'flagstat',
      'base_work_dir'    => $self->o('base_work_dir'),
      'reference_type'   => $self->o('reference_fasta_type'),
      'threads'          => $self->o('samtools_threads'),
      'copy_input'       => $self->o('copy_input_to_temp'),
     },
    -flow_into   => {
        1 => ['samtools_idxstat_summary_for_cellranger'],
      },
  };
  
  ## samtools idxstat metrics
  push @pipeline, {
    -logic_name  => 'samtools_idxstat_summary_for_cellranger',
    -module      => 'ehive.runnable.process.alignment.RunSamtools',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'bam_file'         => '#bam_file#',
      'samtools_command' => 'idxstats',
      'base_work_dir'    => $self->o('base_work_dir'),
      'reference_type'   => $self->o('reference_fasta_type'),
      'copy_input'       => $self->o('copy_input_to_temp'),
     },
    -flow_into   => {
        1 => ['multiqc_report_for_cellranger'],
      },
  };
  
  ## multiqc report building
  push @pipeline, {
    -logic_name  => 'multiqc_report_for_cellranger',
    -module      => 'ehive.runnable.process.alignment.RunAnalysisMultiQC',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'base_results_dir' => $self->o('base_results_dir'),
      'collection_name'  => '#experiment_igf_id#',
      'collection_type'  => $self->o('multiqc_type'),
      'collection_table' => $self->o('cellranger_collection_table'),
      'analysis_name'    => $self->o('multiqc_analysis_name'),
      'tag_name'         => '#species_name#',
      'multiqc_exe'      => $self->o('multiqc_exe'),
      'multiqc_options'  => $self->o('multiqc_options'),
     },
    -flow_into   => {
        1 => ['copy_sample_multiqc_to_remote'],
      },
  };
  
  ## copy multiqc to remote
  push @pipeline, {
      -logic_name   => 'copy_sample_multiqc_to_remote',
      -module       => 'ehive.runnable.process.alignment.CopyAnalysisFilesToRemote',
      -language     => 'python3',
      -meadow_type  => 'PBSPro',
      -rc_name      => '1Gb',
      -analysis_capacity => 2,
      -parameters  => {
        'file_list'           => ['#multiqc_html#'],
        'remote_user'         => $self->o('seqrun_user'),
        'remote_host'         => $self->o('remote_host'),
        'remote_project_path' => $self->o('remote_project_path'),
        },
      -flow_into    => {
          1 => ['mark_experiment_finished'],
      },
  };
  
  ## mark experiment as done
  push @pipeline, {
      -logic_name   => 'mark_experiment_finished',
      -module       => 'ehive.runnable.process.ChangePipelineSeedStatus',
      -language     => 'python3',
      -meadow_type  => 'LOCAL',
      -parameters   => {
        'igf_id'        => '#experiment_igf_id#',
        'task_id'       => '#project_igf_id#',
        'new_status'    => 'FINISHED',
        'pipeline_name' => $self->o('pipeline_name'),
        },
      -flow_into    => {
          1 => ['update_project_analysis'],
      },
  };
  
  ## update analysis page
  push @pipeline, {
      -logic_name   => 'update_project_analysis',
      -module       => 'ehive.runnable.process.UpdateProjectAnalysisStats',
      -language     => 'python3',
      -meadow_type  => 'PBSPro',
      -rc_name      => '1Gb',
      -analysis_capacity => 1,
      -parameters   => {
        'collection_type_list' => [$self->o('multiqc_type'),
                                   $self->o('scanpy_type')],
        'remote_project_path'  => $self->o('remote_project_path'),
        'remote_user'          => $self->o('seqrun_user'),
        'remote_host'          => $self->o('remote_host'),
      },
  };
  
  return \@pipeline;
}

1;