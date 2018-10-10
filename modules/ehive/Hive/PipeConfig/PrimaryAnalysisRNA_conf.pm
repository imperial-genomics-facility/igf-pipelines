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
    'star_analysis_name'   => undef,
    'star_multiqc_type'    => undef,
    'bedGraphToBigWig_path'   => undef,
    'star_collection_table'   => undef,
    'star_genomic_cram_type'  => undef,
    'star_bw_collection_type' => undef,
    ## RSEM
    'rsem_reference_type'  => 'TRANSCRIPTOME_RSEM',
    'rsem_threads'         => 4,
    'rsem_memory_limit'    => 4000,
    'rsem_analysis_name'   => undef,
    'rsem_collection_type' => undef,
    'rsem_collection_table' => undef,
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
        '2->A' => WHEN('#library_source# eq #rna_source#' => ['run_factory_for_rnaseq'],
                       ELSE ['mark_experiment_finished']),
        'A->1' => ['mark_experiment_finished'],
    },
  };
  
  
  ## run factory for genomic and transcriptomic data
  push @pipeline, {
    -logic_name  => 'run_factory_for_rnaseq',
    -module      => 'ehive.runnable.jobfactory.alignment.RunFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -flow_into   => {
        '2->A' => ['fetch_fastq_for_run'],
        'A->1' => ['process_star_bams'],
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
        '1' => ['run_star'],
      },
  };
  
  
  ## collect fastp report
  ## copy report to remote dir
  
  
  ## run star alignment
  push @pipeline, {
    -logic_name  => 'run_star',
    -module      => 'ehive.runnable.process.alignment.RunSTAR',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '42Gb8t',
    -analysis_capacity => 1,
    -parameters  => {
      'star_exe'           => $self->o('star_exe'),
      'r1_read_file'       => '#output_read1#',
      'r2_read_file'       =>  '#output_read2#',
      'output_prefix'      => '#run_igf_id#'.'_'.'#chunk_id#',
      'reference_type'     => $self->o('star_reference_type'),
      'reference_gtf_type' => $self->o('reference_gtf_type'),
      'two_pass_mode'      => $self->o('star_two_pass_mode'),
      'run_thread'         => $self->o('star_run_thread'),
      'star_patameters'    => $self->o('star_patameters'),
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
        1 => ['collect_star_genomic_bam_for_exp','collect_star_transcriptomic_bam_for_exp'],
      },
  };
  
  
  ## collect star genomic bam
  push @pipeline, {
    -logic_name  => 'collect_star_genomic_bam_for_exp',
    -module      => 'ehive.runnable.process.alignment.CollectExpAnalysisChunks',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
       'accu_data' => '#star_aligned_genomic_bam#',
       'output_mode' => 'list',
      },
    -flow_into   => {
        1 => {'picard_merge_and_mark_dup_genomic_bam' => {'star_genomic_bams' => '#exp_chunk_list#'}},
      },
  };
  
  
  ## picard merge and mark duplicate genomic bam
  push @pipeline, {
    -logic_name  => 'picard_merge_and_mark_dup_genomic_bam',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_file'     => '#star_genomic_bams#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'MarkDuplicates',
      'base_work_dir'  => $self->o('base_work_dir'),
      'SORT_ORDER'     => 'coordinate',
     },
    -flow_into => {
          1 => { 'star_genomic_bam_analysis_factory_1' => {'merged_star_genomic_bam' => '#analysis_files[0]#' }},
     },
  };
  
  
  ## star genomic bam analysis factory
  push @pipeline, {
    -logic_name  => 'star_genomic_bam_analysis_factory_1',
    -module      => 'ehive.runnable.jobfactory.alignment.AnalysisFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -parameters  => {
      'file_list' => ['#merged_star_genomic_bam#'],
      },
    -flow_into   => {
        '2->A' => ['convert_star_genomic_bam_to_cram',
                   'star_bigwig'
                  ],
        'A->1' => ['picard_aln_summary_for_star'], 
      },
  };
  
  
  ## convert genomic bam to cram
  push @pipeline, {
    -logic_name  => 'convert_star_genomic_bam_to_cram',
    -module      => 'ehive.runnable.process.alignment.ConvertBamToCram',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
        'bam_file'        => '#bam_file#',
        'base_result_dir' => $self->o('base_results_dir'),
        'threads'         => $self->o('samtools_threads'),
        'samtools_exe'     => $self->o('samtools_exe'),
        'collection_name' => '#experiment_igf_id#',
        'collection_type' => $self->o('star_genomic_cram_type'),
        'collection_table'=> $self->o('star_collection_table'),
        'analysis_name'   => $self->o('star_analysis_name'),
        'tag_name'        => '#species_name#',
        'reference_type'  => $self->o('reference_fasta_type'),
     },
     -flow_into   => {
        1 => ['upload_star_genomic_cram_to_irods'],
      },
  };
  
  
  ## copy star genomic cram to irods
  push @pipeline, {
    -logic_name  => 'upload_star_genomic_cram_to_irods',
    -module      => 'ehive.runnable.process.alignment.UploadAnalysisResultsToIrods',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'file_list'     => '#output_cram_list#',
      'irods_exe_dir' => $self->o('irods_exe_dir'),
      'analysis_name' => $self->o('star_analysis_name'),
      'analysis_dir'  => 'analysis',
      'dir_path_list' => ['#analysis_dir#','#sample_igf_id#','#experiment_igf_id#','#analysis_name#'],
      'file_tag'      => '#sample_igf_id#'.' - '.'#experiment_igf_id#'.' - '.'#analysis_name#'.' - '.'#species_name#',
     },
  };
  
  
  ## star bigwig
  push @pipeline, {
    -logic_name  => 'star_bigwig',
    -module      => 'ehive.runnable.process.alignment.RunSTAR',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '42Gb8t',
    -analysis_capacity => 1,
    -parameters  => {
      'star_exe'           => $self->o('star_exe'),
      input_bam            => '#bam_file#',
      'output_prefix'      => '#run_igf_id#'.'_'.'#chunk_id#',
      'reference_type'     => $self->o('star_reference_type'),
      'reference_gtf_type' => $self->o('reference_gtf_type'),
      'two_pass_mode'      => $self->o('star_two_pass_mode'),
      'run_thread'         => $self->o('star_run_thread'),
      'run_mode'           => 'generate_rna_bigwig',
      'bedGraphToBigWig_path' => $self->o('bedGraphToBigWig_path'),
    },
   -flow_into   => {
        1 => ['load_star_bigwig'],
      },
  };
  
  
  ## Load star bigwig
  push @pipeline, {
    -logic_name  => 'load_star_bigwig',
    -module      => 'ehive.runnable.process.alignment.CollectAnalysisFiles',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'      => '#star_bigwigs#',
      'base_results_dir' => $self->o('base_results_dir'),
      'analysis_name'    => $self->o('star_analysis_name'),
      'collection_name'  => '#experiment_igf_id#',
      'tag_name'         => '#species_name#',
      'collection_type'  => $self->o('star_bw_collection_type'),
      'collection_table' => $self->o('star_collection_table'),
     },
    -flow_into   => {
        1 => ['copy_star_bigwig_to_remote'],
      },
  };
  
  
  ## copy star bigwig to remote
  push @pipeline, {
      -logic_name   => 'copy_star_bigwig_to_remote',
      -module       => 'ehive.runnable.process.alignment.CopyAnalysisFilesToRemote',
      -language     => 'python3',
      -meadow_type  => 'PBSPro',
      -rc_name      => '1Gb',
      -analysis_capacity => 2,
      -parameters  => {
        'file_list'           => '#analysis_output_list#',
        'remote_user'         => $self->o('seqrun_user'),
        'remote_host'         => $self->o('remote_host'),
        'remote_project_path' => $self->o('remote_project_path'),
        },
  };
  
  
  ## picard alignment summary metrics for star
  push @pipeline, {
    -logic_name  => 'picard_aln_summary_for_star',
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
        1 => ['picard_base_dist_summary_for_star'],
      },
  };
  
  
  ## picard base distribution summary metrics
  push @pipeline, {
    -logic_name  => 'picard_base_dist_summary_for_star',
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
        1 => ['picard_gc_bias_summary_for_star'],
      },
  };
  
  ## picard gc bias summary metrics
  push @pipeline, {
    -logic_name  => 'picard_gc_bias_summary_for_star',
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
        1 => ['picard_qual_dist_summary_for_star'],
      },
  };
  
  
  ## picard quality distribution summary metrics
  push @pipeline, {
    -logic_name  => 'picard_qual_dist_summary_for_star',
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
        1 => ['picard_rna_metrics_summary_for_star'],
      },
  };
  
  
  ## picard rna metrics summary metrics
  push @pipeline, {
    -logic_name  => 'picard_rna_metrics_summary_for_star',
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
        1 => ['samtools_flagstat_summary_for_star'],
      },
  };
  
  
  ## samtools flagstat metrics
  push @pipeline, {
    -logic_name  => 'samtools_flagstat_summary_for_star',
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
        1 => ['samtools_idxstat_summary_for_star'],
      },
  };
  
  
  ## samtools idxstat metrics
  push @pipeline, {
    -logic_name  => 'samtools_idxstat_summary_for_star',
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
        1 => ['multiqc_report_for_star'],
      },
  };
  
  
  ## multiqc report building
  push @pipeline, {
    -logic_name  => 'multiqc_report_for_star',
    -module      => 'ehive.runnable.process.alignment.RunAnalysisMultiQC',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'base_results_dir' => $self->o('base_results_dir'),
      'collection_name'  => '#experiment_igf_id#',
      'collection_type'  => $self->o('star_multiqc_type'),
      'collection_table' => $self->o('star_collection_table'),
      'analysis_name'    => $self->o('multiqc_analysis'),
      'tag_name'         => '#species_name#',
      'multiqc_exe'      => $self->o('multiqc_exe'),
      'multiqc_options'  => $self->o('multiqc_options'),
     },
    -flow_into   => {
        1 => ['copy_star_multiqc_to_remote'],
      },
  };
  
  
  ## copy multiqc to remote
  push @pipeline, {
      -logic_name   => 'copy_star_multiqc_to_remote',
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
  };
  
  
  ## collect star transcriptomic bam
  push @pipeline, {
    -logic_name  => 'collect_star_transcriptomic_bam_for_exp',
    -module      => 'ehive.runnable.process.alignment.CollectExpAnalysisChunks',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
       'accu_data' => '#star_aligned_trans_bam#',
       'output_mode' => 'file',
      },
    -flow_into   => {
        1 => {'merge_star_transcriptomic_bams'=>{'star_run_trans_bam_list_file' => '#run_chunk_list_file#'}},
      },
  };
  
  
  ## samtools merge transcriptomic bam
  push @pipeline, {
    -logic_name  => 'merge_star_transcriptomic_bams',
    -module      => 'ehive.runnable.process.alignment.RunSamtools',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'input_file'       => '#star_run_trans_bam_list_file#',
      'samtools_command' => 'merge',
      'samtools_exe'     => $self->o('samtools_exe'),
      'base_work_dir'    => $self->o('base_work_dir'),
      'reference_type'   => $self->o('reference_fasta_type'),
      'threads'          => $self->o('samtools_threads'),
     },
    -flow_into   => {
        1 => {'run_rsem' => {'bam_file' => '#analysis_files#[0]'}},
      },
  };
  
  
  ## run rsem on star transcriptomic bam
  push @pipeline, {
    -logic_name  => 'run_rsem',
    -module      => 'ehive.runnable.process.alignment.RunRSEM',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'input_bam'       => '#bam_file#',
      'samtools_command' => 'merge',
      'base_work_dir'    => $self->o('base_work_dir'),
      'reference_type'   => $self->o('rsem_reference_type'),
      'threads'          => $self->o('rsem_threads'),
      'memory_limit'     => $self->o('rsem_memory_limit'),
     },
    -flow_into   => {
        1 => ['load_rsem_results'],
      },
  };
  
  
  ## Load rsem results
  push @pipeline, {
    -logic_name  => 'load_rsem_results',
    -module      => 'ehive.runnable.process.alignment.CollectAnalysisFiles',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'      => '#rsem_output#',
      'base_results_dir' => $self->o('base_results_dir'),
      'analysis_name'    => $self->o('rsem_analysis_name'),
      'collection_name'  => '#experiment_igf_id#',
      'tag_name'         => '#species_name#',
      'collection_type'  => $self->o('rsem_collection_type'),
      'collection_table' => $self->o('rsem_collection_table'),
     },
    -flow_into   => {
        1 => ['upload_rsem_results_to_irods'],
      },
  };
  
  
  ## copy rsem results to irods
  push @pipeline, {
    -logic_name  => 'upload_rsem_results_to_irods',
    -module      => 'ehive.runnable.process.alignment.UploadAnalysisResultsToIrods',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'file_list'     => '#analysis_output_list#',
      'irods_exe_dir' => $self->o('irods_exe_dir'),
      'analysis_name' => $self->o('rsem_analysis_name'),
      'analysis_dir'  => 'analysis',
      'dir_path_list' => ['#analysis_dir#','#sample_igf_id#','#experiment_igf_id#','#analysis_name#'],
      'file_tag'      => '#sample_igf_id#'.' - '.'#experiment_igf_id#'.' - '.'#analysis_name#'.' - '.'#species_name#',
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
