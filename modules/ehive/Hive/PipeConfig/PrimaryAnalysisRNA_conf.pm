=head1 NAME
    ehive::Hive::PipeConfig::PrimaryAnalysisRNA_conf
=cut


package ehive::Hive::PipeConfig::PrimaryAnalysisRNA_conf;

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
    ## Fetch fastq
    'fastq_collection_type'       => undef,
    'fastq_collection_table'      => undef,
    ## Fastp adapter trimming
    'fastp_exe'            => undef,
    'fastp_options_list'   => ['--qualified_quality_phred 15',
                               '--length_required 15'],
    'fastp_run_thread'     => 4,
    'split_by_lines_count' => 5000000,
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
        2 => WHEN('#library_source# eq #genomic_source# || #library_source# eq #rna_source#' => ['run_factory'],
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
        1 =>  WHEN('#library_source# eq #rna_source# && #fastq_counts# > 0' => ['adapter_trim_without_fastq_split'],
                   ELSE ['mark_experiment_finished']),
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
      'split_fastq' => undef,
      'base_work_dir' => $self->o('base_work_dir'),
      'fastp_exe' => $self->o('fastp_exe'),
      'input_fastq_list' => '#fastq_files#',
    },
    -flow_into   => {
        '2->A' => ['fastq_factory_for_star'],
        'A->1' => ['process_star_bams'],
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
          1 => ['picard_add_rg_tag_to_genomic_bam','picard_add_rg_tag_to_transcriptomic_bam'],
    },
  };
  
  
  ## picard add rg tag to run star genomic bam
  push @pipeline, {
    -logic_name  => 'picard_add_rg_tag_to_genomic_bam',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_file'     => '#star_genomic_bam#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'AddOrReplaceReadGroups',
      'base_work_dir'  => $self->o('base_work_dir'),
      'RGID'           => undef,
      'RGLB'           => undef,
      'RGPL'           => undef,
      'RGPU'           => undef,
      'RGSM'           => undef,
      'RGCN'           => 'Imperial Genomics Facility',
      'SORT_ORDER'     => 'coordinate',
     },
    -flow_into => {
          1 => [ '?accu_name=star_aligned_genomic_bam&accu_address={experiment_igf_id}{seed_date_stamp}[]&accu_input_variable=analysis_files[0]' ],
     },
  };
  
  
  ## picard add rg tag to run star transcriptomic bam
  push @pipeline, {
    -logic_name  => 'picard_add_rg_tag_to_transcriptomic_bam',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_file'     => '#star_transcriptomic_bam#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'AddOrReplaceReadGroups',
      'base_work_dir'  => $self->o('base_work_dir'),
      'RGID'           => undef,
      'RGLB'           => undef,
      'RGPL'           => undef,
      'RGPU'           => undef,
      'RGSM'           => undef,
      'RGCN'           => 'Imperial Genomics Facility',
      'SORT_ORDER'     => 'unsorted',
     },
    -flow_into => {
          1 => [ '?accu_name=star_aligned_trans_bam&accu_address={experiment_igf_id}{seed_date_stamp}[]&accu_input_variable=analysis_files[0]' ],
     },
  };
  
  
  ## process star bams
  push @pipeline, {
    -logic_name  => 'process_star_bams',
    -module      => 'ehive.runnable.IGFBaseProcess',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
       'dataflow_params' => {'finished_star'=>1},
      },
    -flow_into   => {
        1 => ['collect_star_genomic_bam','collect_star_transcriptomic_bam'],
      },
  };
  
  
  ## collect star genomic bam
  push @pipeline, {
    -logic_name  => 'collect_star_genomic_bam',
    -module      => 'ehive.runnable.process.alignment.CollectExpAnalysisChunks',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
       'accu_data' => '#star_aligned_genomic_bam#',
       'output_mode' => 'list',
      },
  };
  
  
  ## collect star transcriptomic bam
  push @pipeline, {
    -logic_name  => 'collect_star_transcriptomic_bam',
    -module      => 'ehive.runnable.process.alignment.CollectExpAnalysisChunks',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
       'accu_data' => '#star_aligned_trans_bam#',
       'output_mode' => 'list',
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
