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
    ## Pipeline
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
    'cleanup_bam_dir'     => 0,
    'cram_type'           => 'ANALYSIS_CRAM',
    'copy_input_to_temp'  => 0,
    'patterned_flow_cell_list' => ['NEXTSEQ','HISEQ4000'],
    ## Irods
    'irods_exe_dir'       => undef,
    ## Java
    'java_exe'            => undef,
    'java_param'          => '-Xmx4g',
    ## Picard
    'picard_jar'          => undef,
    ## MultiQC
    'multiqc_analysis'    => 'multiqc',
    'multiqc_exe'         => undef,
    'multiqc_options'     => '{"--zip-data-dir" : ""}',
    'multiqc_type'        => 'MULTIQC_HTML',
    ## Ref genome
    'reference_fasta_type'=> 'GENOME_FASTA',
    'reference_refFlat'   => 'GENE_REFFLAT',
    'reference_gtf_type'  => 'GENE_GTF',
    ## Cellranger
    'cellranger_exe'              => undef,
    'cellranger_param'            => '{"--nopreflight":"","--disable-ui":"","--jobmode":"pbspro","--localcores":"1","--localmem":"1","--mempercore":"4","--maxjobs":"20"}',
    'cellranger_analysis_name'    => 'cellranger_count',
    'cellranger_collection_table' => 'experiment',
    'cellranger_timeout'          => 43200,
    ## SCANPY
    'scanpy_report_template'      => undef,
    'scanpy_type'                 => 'SCANPY_RESULTS',
    ## Fetch fastq
    'fastq_collection_type'       => undef,
    'fastq_collection_table'      => undef,
    ## Fastp adapter trimming
    'fastp_exe'            => undef,
    'fastp_options_list'   => ['--qualified_quality_phred 15',
                               '--length_required 15'],
    'fastp_run_thread'     => 4,
    'split_by_lines_count' => 5000000,
    ## BWA alignment
    'bwa_exe'              => undef,
    'bwa_reference_type'   => 'GENOME_BWA',
    'bwa_run_thread'       => 4,
    ## Samtools
    'samtools_exe'         => undef,
    'samtools_threads'     => 4,
    ## STAR alignment
    'star_exe'             => undef,
    'star_reference_type'  => 'TRANSCRIPTOME_STAR',
    'star_patameters'      => '{"--outFilterMultimapNmax":"20","--alignSJoverhangMin":"8","--alignSJDBoverhangMin":"1","--outFilterMismatchNmax":"999","--outFilterMismatchNoverReadLmax":"0.04","--alignIntronMin":"20","--alignIntronMax":"1000000,"--alignMatesGapMax":"1000000","--outSAMattributes":"NH HI AS NM MD","--limitBAMsortRAM":"12000000000"}',
    'star_run_thread'      => 8,
    'star_two_pass_mode'   => 1,
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
    -module      => 'ehive.runnable.jobfactory.alignment.RunFactory',
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
                    '#library_source# eq #rna_source# && #fastq_counts# > 0' => ['adapter_trim_without_fastq_split'],
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
      'run_thread' => $self->o('fastp_run_thread'),
      'split_fastq' => 1,
      'base_work_dir' => $self->o('base_work_dir'),
      'fastp_exe' => $self->o('fastp_exe'),
      'input_fastq_list' => '#fastq_files#',
    },
    -flow_into   => {
        2 => ['fastq_factory_for_bwa'],
      },
  };

  ## collect fastp report
  ## copy report to remote dir


  ## fastq factory for bwa
  push @pipeline, {
    -logic_name  => 'fastq_factory_for_bwa',
    -module      => 'ehive.runnable.jobfactory.alignment.FastqAlignmentFactory',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '1Gb',
    -analysis_capacity => 10,
    -parameters  => {
      'read1_list' => '#output_read1#',
      'read2_list' => '#output_read2#',
    },
    -flow_into   => {
        2 => ['run_bwa'],
        
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
      'run_thread' => $self->o('bwa_run_thread'),
       'bwa_exe' => $self->o('bwa_exe'),
       'samtools_exe' => $self->o('samtools_exe'),
       'output_prefix' => '#run_igf_id#'.'_'.'#chunk_id#',
      'reference_type' => $self->o('bwa_reference_type')
    },
    -flow_into => {
          1 => [ '?accu_name=bwa_aligned_bam&accu_address={run_igf_id}{seed_date_stamp}[chunk_id]&accu_input_variable=bwa_bam' ],
        },
  };


  ## adapter trim without fastq splitting
  push @pipeline, {
    -logic_name  => 'adapter_trim_without_fastq_split',
    -module      => 'ehive.runnable.process.alignment.RunFastp',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '1Gb4t',
    -analysis_capacity => 10,
    -parameters  => {
      'fastp_options_list' => $self->o('fastp_options_list'),
      'split_by_lines_count' => $self->o('split_by_lines_count'),
      'run_thread' => $self->o('fastp_run_thread'),
      'base_work_dir' => $self->o('base_work_dir'),
      'fastp_exe' => $self->o('fastp_exe'),
      'input_fastq_list' => '#fastq_files#',
    },
    -flow_into   => {
        2 => ['fastq_factory_for_star'],
      },
  };

  ## collect fastp report
  ## copy report to remote dir


  ## fastq factory for star
  push @pipeline, {
    -logic_name  => 'fastq_factory_for_star',
    -module      => 'ehive.runnable.jobfactory.alignment.FastqAlignmentFactory',
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
      'output_prefix' => '#run_igf_id#'.'_'.'#chunk_id#',
      'reference_type' => $self->o('star_reference_type'),
      'reference_gtf_type' => $self->o('reference_gtf_type'),
      'two_pass_mode' => $self->o('star_two_pass_mode'),
      'run_thread' => $self->o('star_run_thread'),
      'star_patameters' => $self->o('star_patameters'),
    },
    -flow_into => {
          1 => [ '?accu_name=star_aligned_genomic_bam&accu_address={experiment_igf_id}{seed_date_stamp}[]&accu_input_variable=star_genomic_bam' ],
          1 => [ '?accu_name=star_aligned_trans_bam&accu_address={experiment_igf_id}{seed_date_stamp}[]&accu_input_variable=star_transcriptomic_bam' ],
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
  };
  return \@pipeline;
}

1;