=head1 NAME
    ehive::Hive::PipeConfig::PrimaryAnalysisNew_conf
=cut


package ehive::Hive::PipeConfig::PrimaryAnalysisNew_conf;

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
    'genomic_source'      => 'GENOMIC',
    'rna_source'          => 'TRANSCRIPTOMIC',
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
                  '#library_source# eq #genomic_source# || #library_source# eq #rna_source#' => ['run_factory'],
                 ELSE ['mark_experiment_finished'],),
    },
  };

  ## run factory for genomic and transcriptomic data
  push @pipeline, {
    -logic_name  => 'run_factory',
    -module      => 'ehive.runnable.jobfactory.RunFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -flow_into   => {
        2 => ['fetch_fastq_for_run'],
      },
  };

  ## fetch fastq files for a run
  push @pipeline, {
    -logic_name  => 'fetch_fastq_for_run',
    -module      => 'ehive.runnable.process.alignment.FetchFastqForRun',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '1Gb',
    -analysis_capacity => 10,
    -parameters  => {
      'fastq_collection_type' => $self->o('fastq_collection_type'),
      'fastq_collection_table' => $self->o('fastq_collection_table'),
    },
    -flow_into   => {
        1 =>  WHEN('#library_source# eq #genomic_source# && #fastq_counts# > 0 ' => ['adapter_trim_and_fastq_split'],
                    '#library_source# eq #rna_source# && #fastq_counts# > 0' => ['adapter_trim'],
                  ),
      },
  };

  ## adapter trim and fastq splitting by read number
  push @pipeline, {
    -logic_name  => 'adapter_trim_and_fastq_split',
    -module      => 'ehive.runnable.process.alignment.RunFastp',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '1Gb4t',
    -analysis_capacity => 10,
    -parameters  => {
      'fastp_options_list' => $self->o('fastp_options_list'),
      'split_by_lines_count' => $self->o('split_by_lines_count'),
      'run_thread' => 4,
      'split_fastq' => 1,
      'base_work_dir' => $self->o('base_work_dir'),
      'fastp_exe' => $self->o('fastp_exe'),
      'input_fastq_list' => '#fastq_files#',
    },
    -flow_into   => {
        '2->A' => ['fastq_factory_for_bwa'],
        'A->1' => ['merge_run_level_genomic_bams'],
      },
  };

  ## collect fastp report
  ## copy report to remote dir


  ## fastq factory for bwa
  push @pipeline, {
    -logic_name  => 'fastq_factory_for_bwa',
    -module      => 'ehive.runnable.jobfactory.FastqAlignmentFactory',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '1Gb',
    -analysis_capacity => 10,
    -parameters  => {
      'read1_list' => '#output_read1#',
      'read2_list' => '#output_read2#',
    },
    -flow_into   => {
        '2->A' => ['run_bwa'],
        'A->1' => ['merge_chunk_bams'],
      },
  };


  ## run bwa alignment
  push @pipeline, {
    -logic_name  => 'run_bwa',
    -module      => 'ehive.runnable.process.alignment.RunBWA',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 10,
    -parameters  => {
      'run_thread' => 4,
       'bwa_exe' => $self->o('bwa_exe'),
       'samtools_exe' => $self->o('samtools_exe'),
       'output_prefix' => '#run_igf_id#'.'_'.'#chunk_id#',
      'reference_type' => $self->o('bwa_reference_type')
    },
  };


  ## adapter trim without fastq splitting
  push @pipeline, {
    -logic_name  => 'adapter_trim_and_fastq_split',
    -module      => 'ehive.runnable.process.alignment.RunFastp',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '1Gb4t',
    -analysis_capacity => 10,
    -parameters  => {
      'fastp_options_list' => $self->o('fastp_options_list'),
      'split_by_lines_count' => $self->o('split_by_lines_count'),
      'run_thread' => 4,
      'base_work_dir' => $self->o('base_work_dir'),
      'fastp_exe' => $self->o('fastp_exe'),
      'input_fastq_list' => '#fastq_files#',
    },
    -flow_into   => {
        '2->A' => ['fastq_factory_for_star'],
        'A->1' => ['merge_run_level_rna_bams'],
      },
  };

  ## collect fastp report
  ## copy report to remote dir


  ## fastq factory for star
  push @pipeline, {
    -logic_name  => 'fastq_factory_for_bwa',
    -module      => 'ehive.runnable.jobfactory.FastqAlignmentFactory',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '1Gb',
    -analysis_capacity => 10,
    -parameters  => {
      'read1_list' => '#output_read1#',
      'read2_list' => '#output_read2#',
    },
    -flow_into   => {
        2 => ['run_star'],
      },
  };


  ## run star alignment
  push @pipeline, {
    -logic_name  => 'run_star',
    -module      => 'ehive.runnable.process.alignment.RunSTAR',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '42Gb8t',
    -analysis_capacity => 1,
    -parameters  => {
      'star_exe' => $self->o('star_exe'),
      'output_prefix' => '#run_igf_id#',
      'reference_type' => $self->o('star_reference_type'),
      'reference_gtf_type' => $self->o('reference_gtf_type'),
      'two_pass_mode' => 1,
      'run_thread' => 8,
      'star_patameters' => $self->o('star_patameters'),
    },
  };

  ## run cellranger for each experiments
  push @pipeline, {
    -logic_name  => 'run_cellranger_count_for_experiment',
    -module      => 'ehive.runnable.process.alignment.RunCellrangerCount',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'cellranger_exe'     => $self->o('cellranger_exe'),
      'cellranger_options' => $self->o('cellranger_param'),
      'base_work_dir'      => $self->o('base_work_dir'),
      'base_results_dir'   => $self->o('base_results_dir'),
      'job_timeout'        => $self->o('cellranger_timeout'),
      },
  };

  return \@pipeline;
}

1;