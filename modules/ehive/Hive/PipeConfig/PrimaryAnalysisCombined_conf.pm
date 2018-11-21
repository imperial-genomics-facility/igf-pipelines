=head1 NAME
    ehive::Hive::PipeConfig::PrimaryAnalysisCombined_conf
=cut

package ehive::Hive::PipeConfig::PrimaryAnalysisCombined_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('ehive::Hive::PipeConfig::IGFBasePipe_conf');


sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },                                       # here we inherit anything from the base class
    #
    ## PIPELINE
    #---------------------------------------------------------------------------
    'pipeline_name'       => 'PrimaryAnalysisCombined',
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
    'patterned_flow_cell_list'       => ['NEXTSEQ','HISEQ4000'],
    'center_name'         => 'Imperial BRC Genomics Facility',
    #
    ## IRODS
    #---------------------------------------------------------------------------
    'irods_exe_dir'       => undef,
    #
    ## JAVA
    #---------------------------------------------------------------------------
    'java_exe'            => undef,
    'java_param'          => '-Xmx4g',
    #
    ## PICARD
    #---------------------------------------------------------------------------
    'picard_jar'                     => undef,
    'illumina_platform_name'         => 'ILLUMINA',
    #
    ## MULTIQC
    #---------------------------------------------------------------------------
    'multiqc_analysis'    => 'multiqc',
    'multiqc_exe'         => undef,
    'multiqc_options'     => '{"--zip-data-dir" : ""}',
    'multiqc_type'        => 'MULTIQC_HTML',
    'tool_order_list_dnaseq'         => ['fastp','picard','samtools'],
    'tool_order_list_rnaseq'         => ['fastp','star','picard','samtools','featureCounts'],
    'tool_order_list_singlecell'     => ['picard','samtools'],
    'multiqc_template_file'          => undef,
    #
    ## REF GENOME
    #---------------------------------------------------------------------------
    'reference_fasta_type'=> 'GENOME_FASTA',
    'reference_refFlat'   => 'GENE_REFFLAT',
    'reference_gtf_type'  => 'GENE_GTF',
    'two_bit_genome_type' => 'GENOME_TWOBIT_URI',
    #
    ## FASTQ
    #---------------------------------------------------------------------------
    'fastq_collection_type'          => undef,
    'fastq_collection_table'         => undef,
    #
    ## FASTP
    #---------------------------------------------------------------------------
    'fastp_exe'            => undef,
    'fastp_options_list'   => ['--qualified_quality_phred=15','--length_required=15'],
    'fastp_run_thread'     => 4,
    'split_by_lines_count' => 5000000,
    'fastp_analysis_name'  => 'fastp',
    'fastp_html_collection_type'     => 'FASTP_REPORT',
    'fastp_collection_table'         => 'run',
    #
    ## SAMTOOLS
    #---------------------------------------------------------------------------
    'samtools_exe'         => undef,
    'samtools_threads'     => 4,
    #
    ## STAR
    #---------------------------------------------------------------------------
    'star_exe'             => undef,
    'star_reference_type'  => 'TRANSCRIPTOME_STAR',
    'star_patameters'      => '{"--outFilterMultimapNmax":"20","--alignSJoverhangMin":"8","--alignSJDBoverhangMin":"1","--outFilterMismatchNmax":"999","--outFilterMismatchNoverReadLmax":"0.04","--alignIntronMin":"20","--alignIntronMax":"1000000","--alignMatesGapMax":"1000000","--outSAMattributes":"NH HI AS NM MD","--limitBAMsortRAM":"12000000000"}',
    'star_run_thread'      => 8,
    'star_two_pass_mode'   => 1,
    'star_analysis_name'   => undef,
    'bedGraphToBigWig_path'          => undef,
    'star_collection_table'          => undef,
    'star_genomic_cram_type'         => undef,
    'star_bw_collection_type'        => undef,
    #
    ## BWA
    #---------------------------------------------------------------------------
    'bwa_exe'              => undef,
    'bwa_reference_type'   => 'GENOME_BWA',
    'bwa_run_thread'       => 8,
    'bwa_parameters'       => '{"-M":""}',
    'bwa_analysis_name'    => undef,
    'bwa_collection_table' => undef,
    'bwa_genomic_cram_type'          => undef,
    #
    ## RSEM
    #---------------------------------------------------------------------------
    'rsem_exe_dir'         => undef,
    'rsem_reference_type'  => 'TRANSCRIPTOME_RSEM',
    'rsem_analysis_name'   => 'rsem',
    'rsem_threads'         => 8,
    'rsem_memory_limit'    => 4000,
    'rsem_analysis_name'   => undef,
    'rsem_collection_type' => undef,
    'rsem_collection_table'          => undef,
    #
    ## FEATURECOUNTS
    #---------------------------------------------------------------------------
    'featurecounts_exe'    => undef,
    'featurecounts_param'  => undef,
    'featurecounts_threads'          => 4,
    'featurecounts_analysis_name'    => 'featureCounts',
    'featurecounts_collection_type'  => 'FEATURE_COUNTS',
    'featurecounts_collection_table' => 'experiment',
    #
    ## CELLRANGER
    #---------------------------------------------------------------------------
    'cellranger_exe'       => undef,
    'cellranger_param'     => '{"--nopreflight":"","--disable-ui":"","--jobmode":"pbspro","--localcores":"1","--localmem":"1","--mempercore":"4","--maxjobs":"20"}',
    'cellranger_timeout'   => 43200,
    'cellranger_collection_table'    => 'experiment',
    'cellranger_analysis_name'       => 'cellranger_count',
    'cellranger_report_type'         => 'CELLRANGER_REPORT',
    #
    ## SCANPY
    #---------------------------------------------------------------------------
    'scanpy_type'                    => 'SCANPY_RESULTS',
    'scanpy_report_template'         => undef,
    #
    ## DEMULTIPLEXING
    #---------------------------------------------------------------------------
    'demultiplexing_pipeline_name'   => undef,
    #
    ## REMOTE QC PAGE
    #---------------------------------------------------------------------------
    'seqrun_user'          => undef,
    'remote_host'          => undef,
    'remote_project_path'  => undef,
    'analysis_dir'         => 'analysis',
    #
    ## GENOME BROWSER
    #---------------------------------------------------------------------------
    'genome_browser_template_file'   => undef,
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
  
  
  ## GENERIC: collect all experiment seeds
  push @pipeline, {
    -logic_name  => 'find_new_experiment_for_analysis',
    -module      => 'ehive.runnable.jobfactory.PipeseedFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -parameters  => {
      'pipeline_name' => $self->o('pipeline_name'),
      'pipeseed_mode' => $self->o('pipeseed_mode'),
      'rna_source'    => $self->o('rna_source'),
      'dna_source'    => $self->o('genomic_source'),
    },
    -flow_into => {
        2 => WHEN('#library_source# eq #rna_source#' => ['run_factory_for_rnaseq'],
                  '#library_source# eq #dna_source#' => ['run_factory_for_dnaseq'],
                  '#experiment_type# eq #tenx_exp_type# && #library_source# eq #singlecell_source#' => ['run_cellranger_count_for_experiment'],
                       ELSE ['mark_experiment_finished']),
    },
  };
  
  
  #############################  RNA-SEQ START    ##############################
  
  
  ## RNA-SEQ: run factory for genomic and transcriptomic data
  push @pipeline, {
    -logic_name  => 'run_factory_for_rnaseq',
    -module      => 'ehive.runnable.jobfactory.alignment.RunFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -flow_into   => {
        '2->A' => ['fetch_fastq_for_rnaseq_run'],
        'A->1' => ['process_star_bams'],
      },
  };
  
  
  ## RNA-SEQ: fetch fastq files for a run
  push @pipeline, {
    -logic_name  => 'fetch_fastq_for_rnaseq_run',
    -module      => 'ehive.runnable.process.alignment.FetchFastqForRun',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
      'fastq_collection_type'  => $self->o('fastq_collection_type'),
      'fastq_collection_table' => $self->o('fastq_collection_table'),
    },
    -flow_into   => {
        1 => ['adapter_trim_without_fastq_split_rnaseq'],
      },
  };
  
  
  ## RNA-SEQ: adapter trim without fastq splitting
  push @pipeline, {
    -logic_name  => 'adapter_trim_without_fastq_split_rnaseq',
    -module      => 'ehive.runnable.process.alignment.RunFastp',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '1Gb4t',
    -analysis_capacity => 10,
    -parameters  => {
      'fastp_options_list'   => $self->o('fastp_options_list'),
      'split_by_lines_count' => $self->o('split_by_lines_count'),
      'run_thread'           => $self->o('fastp_run_thread'),
      'base_work_dir'        => $self->o('base_work_dir'),
      'fastp_exe'            => $self->o('fastp_exe'),
      'input_fastq_list'     => '#fastq_files_list#',
    },
    -flow_into   => {
        '1' => ['load_fastp_report_rnaseq'],
      },
  };
  
  
  ## RNA-SEQ: collect fastp report
  push @pipeline, {
    -logic_name  => 'load_fastp_report_rnaseq',
    -module      => 'ehive.runnable.process.alignment.CollectAnalysisFiles',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'      => ['#output_html_file#'],
      'base_results_dir' => $self->o('base_results_dir'),
      'analysis_name'    => $self->o('fastp_analysis_name'),
      'collection_name'  => '#run_igf_id#',
      'tag_name'         => '#species_name#',
      'collection_type'  => $self->o('fastp_html_collection_type'),
      'collection_table' => $self->o('fastp_collection_table'),
      'file_suffix'      => 'html',
     },
    -flow_into   => {
        1 => ['run_star',
              '?accu_name=fastp_report_rna&accu_address={experiment_igf_id}{seed_date_stamp}[]&accu_input_variable=output_json_file' ],
      },
  };
  
  
  ## RNA-SEQ: run star alignment
  push @pipeline, {
    -logic_name  => 'run_star',
    -module      => 'ehive.runnable.process.alignment.RunSTAR',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '42Gb8t',
    -analysis_capacity => 10,
    -parameters  => {
      'star_exe'           => $self->o('star_exe'),
      'r1_read_file'       => '#output_read1#',
      'r2_read_file'       => '#output_read2#',
      'output_prefix'      => '#run_igf_id#',
      'base_work_dir'      => $self->o('base_work_dir'),
      'reference_type'     => $self->o('star_reference_type'),
      'reference_gtf_type' => $self->o('reference_gtf_type'),
      'two_pass_mode'      => $self->o('star_two_pass_mode'),
      'run_thread'         => $self->o('star_run_thread'),
      'star_patameters'    => $self->o('star_patameters'),
    },
    -flow_into => {
          1 => ['picard_add_rg_tag_to_star_genomic_bam',
                'picard_add_rg_tag_to_star_transcriptomic_bam',
                '?accu_name=star_logs&accu_address={experiment_igf_id}{seed_date_stamp}[]&accu_input_variable=star_log_file',
                '?accu_name=star_gene_counts&accu_address={experiment_igf_id}{seed_date_stamp}[]&accu_input_variable=star_gene_count_file'
               ]
    },
  };
  
  
  ## RNA-SEQ: picard add rg tag to run star genomic bam
  push @pipeline, {
    -logic_name  => 'picard_add_rg_tag_to_star_genomic_bam',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 10,
    -parameters  => {
      'input_files'    => ['#star_genomic_bam#'],
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'AddOrReplaceReadGroups',
      'base_work_dir'  => $self->o('base_work_dir'),
      'picard_option'  => {
         'RGID'       => '#run_igf_id#',
         'RGLB'       => '#library_name#',
         'RGPL'       => $self->o('illumina_platform_name'),
         'RGPU'       => '#run_igf_id#',
         'RGSM'       => '#sample_igf_id#',
         'RGCN'       => $self->o('center_name'),
         'SORT_ORDER' => 'coordinate',
         },
      'output_prefix'  => '#run_igf_id#'.'_'.'genomic',
     },
    -flow_into => {
          1 => [ '?accu_name=star_aligned_genomic_bam&accu_address={experiment_igf_id}{seed_date_stamp}[]&accu_input_variable=analysis_files' ],
     },
  };
  
  
  ## RNA-SEQ: picard add rg tag to run star transcriptomic bam
  push @pipeline, {
    -logic_name  => 'picard_add_rg_tag_to_star_transcriptomic_bam',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 10,
    -parameters  => {
      'input_files'    => ['#star_transcriptomic_bam#'],
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'AddOrReplaceReadGroups',
      'base_work_dir'  => $self->o('base_work_dir'),
      'picard_option'  => {
         'RGID'       => '#run_igf_id#',
         'RGLB'       => '#library_name#',
         'RGPL'       => $self->o('illumina_platform_name'),
         'RGPU'       => '#run_igf_id#',
         'RGSM'       => '#sample_igf_id#',
         'RGCN'       => $self->o('center_name'),
         'SORT_ORDER' => 'unsorted',
         },
      'output_prefix'  => '#run_igf_id#'.'_'.'transcriptomic',
     },
    -flow_into => {
          1 => [ '?accu_name=star_aligned_trans_bam&accu_address={experiment_igf_id}{seed_date_stamp}[]&accu_input_variable=analysis_files' ],
     },
  };
  
  
  ## RNA-SEQ: process star bams
  push @pipeline, {
    -logic_name  => 'process_star_bams',
    -module      => 'ehive.runnable.IGFBaseJobFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
       'sub_tasks' => [{'pseudo_exp_id'=> '#experiment_igf_id#'}],
      },
    -flow_into   => {
        '2->A' => ['collect_star_genomic_bam_for_exp',
                   'collect_star_transcriptomic_bam_for_exp'],
        'A->1' => ['mark_experiment_finished'],
      },
  };
  
  
  ## RNA-SEQ: collect star genomic bam
  push @pipeline, {
    -logic_name  => 'collect_star_genomic_bam_for_exp',
    -module      => 'ehive.runnable.process.alignment.CollectExpAnalysisChunks',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
       'accu_data'     => '#star_aligned_genomic_bam#',
       'output_mode'   => 'list',
       'base_work_dir' => $self->o('base_work_dir'),
      },
    -flow_into   => {
        1 => {'picard_merge_and_mark_dup_star_genomic_bam' => {'star_genomic_bams' => '#exp_chunk_list#'}},
      },
  };
  
  
  ## RNA-SEQ: picard merge and mark duplicate genomic bam
  push @pipeline, {
    -logic_name  => 'picard_merge_and_mark_dup_star_genomic_bam',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '8Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => '#star_genomic_bams#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'MarkDuplicates',
      'base_work_dir'  => $self->o('base_work_dir'),
      'picard_option'  => { 'ASSUME_SORT_ORDER' => 'coordinate'},
      'output_prefix'  => '#experiment_igf_id#',
     },
    -flow_into => {
          1 => { 'star_genomic_bam_analysis_factory' => {'merged_star_genomic_bams' => '#bam_files#',
                                                         'analysis_files' => '#analysis_files#'}},
     },
  };
  
  
  ## RNA-SEQ: star genomic bam analysis factory
  push @pipeline, {
    -logic_name  => 'star_genomic_bam_analysis_factory',
    -module      => 'ehive.runnable.jobfactory.alignment.AnalysisFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -parameters  => {
      'file_list' => '#merged_star_genomic_bams#',
      },
    -flow_into   => {
        '2->A' => ['convert_star_genomic_bam_to_cram',
                   'star_bigwig'
                  ],
        'A->1' => {'collect_star_log_for_exp' => {'merged_star_genomic_bams' => '#merged_star_genomic_bams#',
                                                  'analysis_files' => '#analysis_files#'}},
      },
  };
  
  
  ## RNA-SEQ: convert genomic bam to cram
  push @pipeline, {
    -logic_name  => 'convert_star_genomic_bam_to_cram',
    -module      => 'ehive.runnable.process.alignment.ConvertBamToCram',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
        'bam_files'       => ['#input_file#'],                                  # updated
        'base_result_dir' => $self->o('base_results_dir'),
        'threads'         => $self->o('samtools_threads'),
        'samtools_exe'    => $self->o('samtools_exe'),
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
  
  
  ## RNA-SEQ: copy star genomic cram to irods
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
     -flow_into   => {
        1 => ['run_featureCounts_for_rnaseq'],
      },
  };
  
  
  ## RNA-SEQ: run feature counts
  push @pipeline, {
    -logic_name  => 'run_featureCounts_for_rnaseq',
    -module      => 'ehive.runnable.process.alignment.RunFeatureCounts',
    -language    => 'python3',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'featurecounts_exe'  => $self->o('featurecounts_exe'),
      'input_files'        => ['#input_file#'],
      'reference_gtf'      => $self->o('reference_gtf_type'),
      'run_thread'         => $self->o('featurecounts_threads'),
      'parameter_options'  => $self->o('featurecounts_param'),
      'output_prefix'      => '#experiment_igf_id#',
      'base_work_dir'      => $self->o('base_work_dir'),
    },
    -flow_into   => {
        1 => ['load_featurecounts_results',
              '?accu_name=feature_count_logs&accu_address={experiment_igf_id}{seed_date_stamp}[]&accu_input_variable=featureCounts_summary'
             ],
      },
  };
  
  
  ## RNA-SEQ: load featureCounts results
  push @pipeline, {
    -logic_name  => 'load_featurecounts_results',
    -module      => 'ehive.runnable.process.alignment.CollectAnalysisFiles',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'      => ['#featureCounts_output#'],
      'base_results_dir' => $self->o('base_results_dir'),
      'analysis_name'    => $self->o('featurecounts_analysis_name'),
      'collection_name'  => '#experiment_igf_id#',
      'tag_name'         => '#species_name#',
      'file_suffix'      => 'txt',
      'collection_type'  => $self->o('featurecounts_collection_type'),
      'collection_table' => $self->o('featurecounts_collection_table'),
     },
    -flow_into   => {
        1 => ['upload_featurecounts_results_to_irods'],
      },
  };
  
  
  ## RNA-SEQ: copy featureCounts results to irods
  push @pipeline, {
    -logic_name  => 'upload_featurecounts_results_to_irods',
    -module      => 'ehive.runnable.process.alignment.UploadAnalysisResultsToIrods',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'file_list'     => '#analysis_output_list#',
      'irods_exe_dir' => $self->o('irods_exe_dir'),
      'analysis_name' => $self->o('featurecounts_analysis_name'),
      'analysis_dir'  => 'analysis',
      'dir_path_list' => ['#analysis_dir#','#sample_igf_id#','#experiment_igf_id#','#analysis_name#'],
      'file_tag'      => '#sample_igf_id#'.' - '.'#experiment_igf_id#'.' - '.'#analysis_name#'.' - '.'#species_name#',
     },
  };
  
  
  ## RNA-SEQ: star bigwig
  push @pipeline, {
    -logic_name  => 'star_bigwig',
    -module      => 'ehive.runnable.process.alignment.RunSTAR',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '42Gb8t',
    -analysis_capacity => 1,
    -parameters  => {
      'star_exe'           => $self->o('star_exe'),
      'input_bam'          => '#input_file#',
      'output_prefix'      => '#experiment_igf_id#',
      'base_work_dir'      => $self->o('base_work_dir'),
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
  
  
  ## RNA-SEQ: load star bigwig
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
  
  
  ## RNA-SEQ: copy star bigwig to remote
  push @pipeline, {
      -logic_name   => 'copy_star_bigwig_to_remote',
      -module       => 'ehive.runnable.process.alignment.CopyAnalysisFilesToRemote',
      -language     => 'python3',
      -meadow_type  => 'PBSPro',
      -rc_name      => '1Gb',
      -analysis_capacity => 2,
      -parameters  => {
        'analysis_dir'        => $self->o('analysis_dir'),
        'dir_labels'          => ['#analysis_dir#','#sample_igf_id#'],
        'file_list'           => '#analysis_output_list#',
        'remote_user'         => $self->o('seqrun_user'),
        'remote_host'         => $self->o('remote_host'),
        'remote_project_path' => $self->o('remote_project_path'),
        },
  };
  
  
  ## RNA-SEQ: collect star logs
  push @pipeline, {
    -logic_name  => 'collect_star_log_for_exp',
    -module      => 'ehive.runnable.process.alignment.CollectExpAnalysisChunks',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
       'exp_chunk_list'  => '#analysis_files#',
       'accu_data'       => '#star_logs#',
       'output_mode'     => 'list',
       'base_work_dir'   => $self->o('base_work_dir'),
      },
    -flow_into   => {
        1 => ['collect_fastp_json_for_exp_rnaseq'],
      },
  };
  
  
  ## RNA-SEQ: collect fastp json
  push @pipeline, {
    -logic_name  => 'collect_fastp_json_for_exp_rnaseq',
    -module      => 'ehive.runnable.process.alignment.CollectExpAnalysisChunks',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
       'accu_data'     => '#fastp_report_rna#',
       'output_mode'   => 'list',
       'base_work_dir' => $self->o('base_work_dir'),
      },
    -flow_into   => {
        1 => ['collect_featureCounts_for_exp_rnaseq'],
      },
  };
  
  
  ## RNA-SEQ: collect featureCounts summary
  push @pipeline, {
    -logic_name  => 'collect_featureCounts_for_exp_rnaseq',
    -module      => 'ehive.runnable.process.alignment.CollectExpAnalysisChunks',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
       'accu_data'     => '#feature_count_logs#',
       'output_mode'   => 'list',
       'base_work_dir' => $self->o('base_work_dir'),
      },
    -flow_into   => {
        1 => {'picard_aln_summary_for_star' => {'analysis_files' => '#exp_chunk_list#'}},
      },
  };
  
  
  ## RNA-SEQ: picard alignment summary metrics for star
  push @pipeline, {
    -logic_name  => 'picard_aln_summary_for_star',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => '#merged_star_genomic_bams#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'CollectAlignmentSummaryMetrics',
      'base_work_dir'  => $self->o('base_work_dir'),
      'reference_type' => $self->o('reference_fasta_type'),
      'output_prefix'  => '#experiment_igf_id#',
     },
    -flow_into   => {
        1 => ['picard_base_dist_summary_for_star'],
      },
  };
  
  
  ## RNA-SEQ: picard base distribution summary metrics
  push @pipeline, {
    -logic_name  => 'picard_base_dist_summary_for_star',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => '#merged_star_genomic_bams#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'CollectBaseDistributionByCycle',
      'base_work_dir'  => $self->o('base_work_dir'),
      'reference_type' => $self->o('reference_fasta_type'),
      'output_prefix'  => '#experiment_igf_id#',
     },
    -flow_into   => {
        1 => ['picard_gc_bias_summary_for_star'],
      },
  };
  
  
  ## RNA-SEQ: picard gc bias summary metrics
  push @pipeline, {
    -logic_name  => 'picard_gc_bias_summary_for_star',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => '#merged_star_genomic_bams#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'CollectGcBiasMetrics',
      'base_work_dir'  => $self->o('base_work_dir'),
      'reference_type' => $self->o('reference_fasta_type'),
      'output_prefix'  => '#experiment_igf_id#',
     },
    -flow_into   => {
        1 => ['picard_qual_dist_summary_for_star'],
      },
  };
  
  
  ## RNA-SEQ: picard quality distribution summary metrics
  push @pipeline, {
    -logic_name  => 'picard_qual_dist_summary_for_star',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => '#merged_star_genomic_bams#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'QualityScoreDistribution',
      'base_work_dir'  => $self->o('base_work_dir'),
      'reference_type' => $self->o('reference_fasta_type'),
      'output_prefix'  => '#experiment_igf_id#',
     },
    -flow_into   => {
        1 => ['picard_rna_metrics_summary_for_star'],
      },
  };
  
  
  ## RNA-SEQ: picard rna metrics summary metrics
  push @pipeline, {
    -logic_name  => 'picard_rna_metrics_summary_for_star',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => '#merged_star_genomic_bams#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'CollectRnaSeqMetrics',
      'base_work_dir'  => $self->o('base_work_dir'),
      'reference_type' => $self->o('reference_fasta_type'),
      'reference_refFlat' => $self->o('reference_refFlat'),
      'output_prefix'  => '#experiment_igf_id#',
     },
    -flow_into   => {
        1 => ['samtools_stats_summary_for_star'],
      },
  };
  
  
  ## RNA-SEQ: samtools stats metrics
  push @pipeline, {
    -logic_name  => 'samtools_stats_summary_for_star',
    -module      => 'ehive.runnable.process.alignment.RunSamtools',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'      => '#merged_star_genomic_bams#',
      'samtools_command' => 'stats',
      'output_prefix'    => '#experiment_igf_id#',
      'base_work_dir'    => $self->o('base_work_dir'),
      'reference_type'   => $self->o('reference_fasta_type'),
      'samtools_exe'     => $self->o('samtools_exe'),
      'threads'          => $self->o('samtools_threads'),
     },
    -flow_into   => {
        1 => ['samtools_idxstat_summary_for_star'],
      },
  };
  
  
  ## RNA-SEQ: samtools idxstat metrics
  push @pipeline, {
    -logic_name  => 'samtools_idxstat_summary_for_star',
    -module      => 'ehive.runnable.process.alignment.RunSamtools',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'      => '#merged_star_genomic_bams#',
      'samtools_command' => 'idxstats',
      'output_prefix'    => '#experiment_igf_id#',
      'base_work_dir'    => $self->o('base_work_dir'),
      'samtools_exe'     => $self->o('samtools_exe'),
      'reference_type'   => $self->o('reference_fasta_type'),
     },
    -flow_into   => {
        1 => ['multiqc_report_for_star'],
      },
  };
  
  
  ## RNA-SEQ: multiqc report building
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
      'tool_order_list'  => $self->o('tool_order_list_rnaseq'),
      'multiqc_template_file'  => $self->o('multiqc_template_file'),
     },
    -flow_into   => {
        1 => ['copy_star_multiqc_to_remote'],
      },
  };
  
  
  ## RNA-SEQ: copy multiqc to remote
  push @pipeline, {
      -logic_name   => 'copy_star_multiqc_to_remote',
      -module       => 'ehive.runnable.process.alignment.CopyAnalysisFilesToRemote',
      -language     => 'python3',
      -meadow_type  => 'PBSPro',
      -rc_name      => '1Gb',
      -analysis_capacity => 2,
      -parameters  => {
        'analysis_dir'        => $self->o('analysis_dir'),
        'dir_labels'          => ['#analysis_dir#','#sample_igf_id#'],
        'file_list'           => ['#multiqc_html#'],
        'remote_user'         => $self->o('seqrun_user'),
        'remote_host'         => $self->o('remote_host'),
        'remote_project_path' => $self->o('remote_project_path'),
        },
  };
  
  
  ## RNA-SEQ: collect star transcriptomic bam
  push @pipeline, {
    -logic_name  => 'collect_star_transcriptomic_bam_for_exp',
    -module      => 'ehive.runnable.process.alignment.CollectExpAnalysisChunks',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
       'accu_data'      => '#star_aligned_trans_bam#',
       'output_mode'    => 'file',
       'base_work_dir'  => $self->o('base_work_dir'),
      },
    -flow_into   => {
        1 => {'merge_star_transcriptomic_bams'=>{'star_run_trans_bam_list_file' => '#exp_chunk_list_file#'}},
      },
  };
  
  
  ## RNA-SEQ: samtools merge transcriptomic bam
  push @pipeline, {
    -logic_name  => 'merge_star_transcriptomic_bams',
    -module      => 'ehive.runnable.process.alignment.RunSamtools',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'      => ['#star_run_trans_bam_list_file#'],
      'samtools_command' => 'merge',
      'output_prefix'    => '#experiment_igf_id#',
      'sorted_by_name'   => 1,                                                  # enable read sorting by name for rsem
      'samtools_exe'     => $self->o('samtools_exe'),
      'base_work_dir'    => $self->o('base_work_dir'),
      'reference_type'   => $self->o('reference_fasta_type'),
      'threads'          => $self->o('samtools_threads'),
     },
    -flow_into   => {
        1 => {'run_rsem' => {'bam_files' => '#analysis_files#'}},
      },
  };
  
  
  ## RNA-SEQ: run rsem on star transcriptomic bam
  push @pipeline, {
    -logic_name  => 'run_rsem',
    -module      => 'ehive.runnable.process.alignment.RunRSEM',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '16Gb8t',
    -analysis_capacity => 2,
    -parameters  => {
      'input_bams'     => '#bam_files#',
      'rsem_exe_dir'   => $self->o('rsem_exe_dir'),
      'base_work_dir'  => $self->o('base_work_dir'),
      'output_prefix'  => '#experiment_igf_id#',
      'species_name'   => '#species_name#',
      'reference_type' => $self->o('rsem_reference_type'),
      'threads'        => $self->o('rsem_threads'),
      'memory_limit'   => $self->o('rsem_memory_limit'),
     },
    -flow_into   => {
        1 => ['load_rsem_results',
              '?accu_name=rsem_logs&accu_address={experiment_igf_id}{seed_date_stamp}[]&accu_input_variable=rsem_log_file' ],
      },
  };
  
  
  ## RNA-SEQ: load rsem results
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
  
  
  ## RNA-SEQ: copy rsem results to irods
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
      'analysis_dir'  => $self->o('analysis_dir'),
      'dir_path_list' => ['#analysis_dir#','#sample_igf_id#','#experiment_igf_id#','#analysis_name#'],
      'file_tag'      => '#sample_igf_id#'.' - '.'#experiment_igf_id#'.' - '.'#analysis_name#'.' - '.'#species_name#',
     },
  };
  
  
  #############################  RNA-SEQ END      ##############################
  #############################  DNA-SEQ START    ##############################
  
  
  ## DNA-SEQ: run factory for genomic data
  push @pipeline, {
    -logic_name  => 'run_factory_for_dnaseq',
    -module      => 'ehive.runnable.jobfactory.alignment.RunFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -flow_into   => {
        '2->A' => ['fetch_fastq_for_dnaseq_run'],
        'A->1' => ['process_bwa_bams'],
      },
  };
  
  
  ## DNA-SEQ: fetch fastq files for a run
  push @pipeline, {
    -logic_name  => 'fetch_fastq_for_dnaseq_run',
    -module      => 'ehive.runnable.process.alignment.FetchFastqForRun',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
      'fastq_collection_type'  => $self->o('fastq_collection_type'),
      'fastq_collection_table' => $self->o('fastq_collection_table'),
    },
    -flow_into   => {
        1 => ['adapter_trim_without_fastq_split_dnaseq'],
      },
  };
  
  
  ## DNA-SEQ: adapter trim without fastq splitting
  push @pipeline, {
    -logic_name  => 'adapter_trim_without_fastq_split_dnaseq',
    -module      => 'ehive.runnable.process.alignment.RunFastp',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '1Gb4t',
    -analysis_capacity => 10,
    -parameters  => {
      'fastp_options_list'   => $self->o('fastp_options_list'),
      'split_by_lines_count' => $self->o('split_by_lines_count'),
      'run_thread'           => $self->o('fastp_run_thread'),
      'base_work_dir'        => $self->o('base_work_dir'),
      'fastp_exe'            => $self->o('fastp_exe'),
      'input_fastq_list'     => '#fastq_files_list#',
    },
    -flow_into   => {
        '1' => ['load_fastp_report_dnaseq'],
      },
  };
  
  
  ## DNA-SEQ: collect fastp report
  push @pipeline, {
    -logic_name  => 'load_fastp_report_dnaseq',
    -module      => 'ehive.runnable.process.alignment.CollectAnalysisFiles',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'      => ['#output_html_file#'],
      'base_results_dir' => $self->o('base_results_dir'),
      'analysis_name'    => $self->o('fastp_analysis_name'),
      'collection_name'  => '#run_igf_id#',
      'tag_name'         => '#species_name#',
      'collection_type'  => $self->o('fastp_html_collection_type'),
      'collection_table' => $self->o('fastp_collection_table'),
      'file_suffix'      => 'html',
     },
    -flow_into   => {
        1 => ['run_bwa',
              '?accu_name=bwa_fastp_report&accu_address={experiment_igf_id}{seed_date_stamp}[]&accu_input_variable=output_json_file'],
      },
  };
  
  
  ## DNA-SEQ: run bwa alignment
  push @pipeline, {
    -logic_name  => 'run_bwa',
    -module      => 'ehive.runnable.process.alignment.RunBWA',                  # FIX ME
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb8t',
    -analysis_capacity => 10,
    -parameters  => {
      'bwa_exe'            => $self->o('bwa_exe'),
      'samtools_exe'       => $self->o('samtools_exe'),
      'r1_read_file'       => '#output_read1#',
      'r2_read_file'       => '#output_read2#',
      'output_prefix'      => '#run_igf_id#',
      'base_work_dir'      => $self->o('base_work_dir'),
      'reference_type'     => $self->o('bwa_reference_type'),
      'run_thread'         => $self->o('bwa_run_thread'),
      'parameter_options'  => $self->o('bwa_parameters'),
    },
    -flow_into => {
          1 => ['picard_add_rg_tag_to_bwa_genomic_bam'],
    },
  };
  
  
  ## DNA-SEQ: picard add rg tag to run bwa genomic bam
  push @pipeline, {
    -logic_name  => 'picard_add_rg_tag_to_bwa_genomic_bam',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 10,
    -parameters  => {
      'input_files'    => ['#bwa_bam#'],
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'AddOrReplaceReadGroups',
      'base_work_dir'  => $self->o('base_work_dir'),
      'picard_option'  => {
         'RGID'       => '#run_igf_id#',
         'RGLB'       => '#library_name#',
         'RGPL'       => $self->o('illumina_platform_name'),
         'RGPU'       => '#run_igf_id#',
         'RGSM'       => '#sample_igf_id#',
         'RGCN'       => $self->o('center_name'),
         'SORT_ORDER' => 'coordinate',
         },
      'output_prefix'  => '#run_igf_id#'.'_'.'genomic',
     },
    -flow_into => {
          1 => [ '?accu_name=bwa_aligned_genomic_bam&accu_address={experiment_igf_id}{seed_date_stamp}[]&accu_input_variable=analysis_files' ],
     },
  };
  
  
  ## DNA-SEQ: process bwa bams
  push @pipeline, {
    -logic_name  => 'process_bwa_bams',
    -module      => 'ehive.runnable.IGFBaseJobFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
       'sub_tasks' => [{'pseudo_exp_id'=> '#experiment_igf_id#'}],
      },
    -flow_into   => {
        '2->A' => ['collect_bwa_genomic_bam_for_exp'],
        'A->1' => ['mark_experiment_finished'],
      },
  };
  
  
  ## DNA-SEQ: collect bwa genomic bam
  push @pipeline, {
    -logic_name  => 'collect_bwa_genomic_bam_for_exp',
    -module      => 'ehive.runnable.process.alignment.CollectExpAnalysisChunks',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
       'accu_data'     => '#bwa_aligned_genomic_bam#',
       'output_mode'   => 'list',
       'base_work_dir' => $self->o('base_work_dir'),
      },
    -flow_into   => {
        1 => {'collect_bwa_fastp_json_for_exp' => {'bwa_genomic_bams' => '#exp_chunk_list#'}},
      },
  };
  
  
  ## DNA-SEQ: collect fastp json
  push @pipeline, {
    -logic_name  => 'collect_bwa_fastp_json_for_exp',
    -module      => 'ehive.runnable.process.alignment.CollectExpAnalysisChunks',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -parameters  => {
       'accu_data'     => '#bwa_fastp_report#',
       'output_mode'   => 'list',
       'base_work_dir' => $self->o('base_work_dir'),
      },
    -flow_into   => {
        1 => {'picard_bwa_merge_and_mark_dup_bwa_genomic_bam' => {'analysis_files' => '#exp_chunk_list#'}},
      },
  };
  
  
  ## DNA-SEQ: picard merge and mark duplicate genomic bam
  push @pipeline, {
    -logic_name  => 'picard_bwa_merge_and_mark_dup_bwa_genomic_bam',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '8Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => '#bwa_genomic_bams#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'MarkDuplicates',
      'base_work_dir'  => $self->o('base_work_dir'),
      'picard_option'  => { 'ASSUME_SORT_ORDER' => 'coordinate'},
      'output_prefix'  => '#experiment_igf_id#',
     },
    -flow_into => {
          1 => { 'convert_bwa_genomic_bam_to_cram' => {'merged_bwa_genomic_bams' => '#bam_files#',
                                                       'analysis_files' => '#analysis_files#'}},
     },
  };
  
  
  ## DNA-SEQ: convert genomic bam to cram
  push @pipeline, {
    -logic_name  => 'convert_bwa_genomic_bam_to_cram',
    -module      => 'ehive.runnable.process.alignment.ConvertBamToCram',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
        'bam_files'       => '#merged_bwa_genomic_bams#',                       # fixed now
        'base_result_dir' => $self->o('base_results_dir'),
        'threads'         => $self->o('samtools_threads'),
        'samtools_exe'    => $self->o('samtools_exe'),
        'collection_name' => '#experiment_igf_id#',
        'collection_type' => $self->o('bwa_genomic_cram_type'),
        'collection_table'=> $self->o('bwa_collection_table'),
        'analysis_name'   => $self->o('bwa_analysis_name'),
        'tag_name'        => '#species_name#',
        'reference_type'  => $self->o('reference_fasta_type'),
     },
     -flow_into   => {
        1 => ['upload_bwa_genomic_cram_to_irods'],
      },
  };
  
  
  ## DNA-SEQ: copy bwa genomic cram to irods
  push @pipeline, {
    -logic_name  => 'upload_bwa_genomic_cram_to_irods',
    -module      => 'ehive.runnable.process.alignment.UploadAnalysisResultsToIrods',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'file_list'     => '#output_cram_list#',
      'irods_exe_dir' => $self->o('irods_exe_dir'),
      'analysis_name' => $self->o('bwa_analysis_name'),
      'analysis_dir'  => 'analysis',
      'dir_path_list' => ['#analysis_dir#','#sample_igf_id#','#experiment_igf_id#','#analysis_name#'],
      'file_tag'      => '#sample_igf_id#'.' - '.'#experiment_igf_id#'.' - '.'#analysis_name#'.' - '.'#species_name#',
     },
     -flow_into   => {
        1 => ['picard_aln_summary_for_bwa'],
      },
  };
  
  
  ## DNA-SEQ: picard alignment summary metrics for bwa
  push @pipeline, {
    -logic_name  => 'picard_aln_summary_for_bwa',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => '#merged_bwa_genomic_bams#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'CollectAlignmentSummaryMetrics',
      'base_work_dir'  => $self->o('base_work_dir'),
      'reference_type' => $self->o('reference_fasta_type'),
      'output_prefix'  => '#experiment_igf_id#',
     },
    -flow_into   => {
        1 => ['picard_base_dist_summary_for_bwa'],
      },
  };
  
  
  ## DNA-SEQ: picard base distribution summary metrics
  push @pipeline, {
    -logic_name  => 'picard_base_dist_summary_for_bwa',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => '#merged_bwa_genomic_bams#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'CollectBaseDistributionByCycle',
      'base_work_dir'  => $self->o('base_work_dir'),
      'reference_type' => $self->o('reference_fasta_type'),
      'output_prefix'  => '#experiment_igf_id#',
     },
    -flow_into   => {
        1 => ['picard_gc_bias_summary_for_bwa'],
      },
  };
  
  
  ## DNA-SEQ: picard gc bias summary metrics
  push @pipeline, {
    -logic_name  => 'picard_gc_bias_summary_for_bwa',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => '#merged_bwa_genomic_bams#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'CollectGcBiasMetrics',
      'base_work_dir'  => $self->o('base_work_dir'),
      'reference_type' => $self->o('reference_fasta_type'),
      'output_prefix'  => '#experiment_igf_id#',
     },
    -flow_into   => {
        1 => ['picard_qual_dist_summary_for_bwa'],
      },
  };
  
  
  ## DNA-SEQ: picard quality distribution summary metrics
  push @pipeline, {
    -logic_name  => 'picard_qual_dist_summary_for_bwa',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => '#merged_bwa_genomic_bams#',
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'QualityScoreDistribution',
      'base_work_dir'  => $self->o('base_work_dir'),
      'reference_type' => $self->o('reference_fasta_type'),
      'output_prefix'  => '#experiment_igf_id#',
     },
    -flow_into   => {
        1 => ['samtools_stats_summary_for_bwa'],
      },
  };
  
  
  ## DNA-SEQ: samtools stats metrics
  push @pipeline, {
    -logic_name  => 'samtools_stats_summary_for_bwa',
    -module      => 'ehive.runnable.process.alignment.RunSamtools',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'      => '#merged_bwa_genomic_bams#',
      'samtools_command' => 'stats',
      'output_prefix'    => '#experiment_igf_id#',
      'base_work_dir'    => $self->o('base_work_dir'),
      'reference_type'   => $self->o('reference_fasta_type'),
      'samtools_exe'     => $self->o('samtools_exe'),
      'threads'          => $self->o('samtools_threads'),
     },
    -flow_into   => {
        1 => ['samtools_idxstat_summary_for_bwa'],
      },
  };
  
  
  ## DNA-SEQ: samtools idxstat metrics
  push @pipeline, {
    -logic_name  => 'samtools_idxstat_summary_for_bwa',
    -module      => 'ehive.runnable.process.alignment.RunSamtools',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'      => '#merged_bwa_genomic_bams#',
      'samtools_command' => 'idxstats',
      'output_prefix'    => '#experiment_igf_id#',
      'base_work_dir'    => $self->o('base_work_dir'),
      'samtools_exe'     => $self->o('samtools_exe'),
      'reference_type'   => $self->o('reference_fasta_type'),
     },
    -flow_into   => {
        1 => ['multiqc_report_for_bwa'],
      },
  };
  
  
  ## DNA-SEQ: multiqc report building
  push @pipeline, {
    -logic_name  => 'multiqc_report_for_bwa',
    -module      => 'ehive.runnable.process.alignment.RunAnalysisMultiQC',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'base_results_dir' => $self->o('base_results_dir'),
      'collection_name'  => '#experiment_igf_id#',
      'collection_type'  => $self->o('multiqc_type'),
      'collection_table' => $self->o('bwa_collection_table'),
      'analysis_name'    => $self->o('multiqc_analysis'),
      'tag_name'         => '#species_name#',
      'multiqc_exe'      => $self->o('multiqc_exe'),
      'multiqc_options'  => $self->o('multiqc_options'),
      'tool_order_list'  => $self->o('tool_order_list_dnaseq'),
      'multiqc_template_file'  => $self->o('multiqc_template_file'),
     },
    -flow_into   => {
        1 => ['copy_bwa_multiqc_to_remote'],
      },
  };
  
  
  ## DNA-SEQ: copy multiqc to remote
  push @pipeline, {
      -logic_name   => 'copy_bwa_multiqc_to_remote',
      -module       => 'ehive.runnable.process.alignment.CopyAnalysisFilesToRemote',
      -language     => 'python3',
      -meadow_type  => 'PBSPro',
      -rc_name      => '1Gb',
      -analysis_capacity => 2,
      -parameters  => {
        'analysis_dir'        => $self->o('analysis_dir'),
        'dir_labels'          => ['#analysis_dir#','#sample_igf_id#'],
        'file_list'           => ['#multiqc_html#'],
        'remote_user'         => $self->o('seqrun_user'),
        'remote_host'         => $self->o('remote_host'),
        'remote_project_path' => $self->o('remote_project_path'),
        },
  };
  
  
  #############################  DNA-SEQ END      ##############################
  #############################  SINGLECELL START ##############################
  
  
  ## SINGLECELL: run cellranger for each experiments
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
  
  
  ## SINGLECELL: load cellranger output for each experiments
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
  
  
  ## SINGLECELL: single cell analysis factory 1
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
                   'convert_cellranger_bam_to_cram',
                   'scanpy_report_generation'
                   ],
        'A->1' => ['picard_aln_summary_for_cellranger'], 
      },
  };
  
  
  ## SINGLECELL: upload cellranger resilts to irods
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
      'analysis_dir'  => $self->o('analysis_dir'),
      'dir_path_list' => ['#analysis_dir#','#sample_igf_id#','#experiment_igf_id#','#analysis_name#'],
      'file_tag'      => '#sample_igf_id#'.' - '.'#experiment_igf_id#'.' - '.'#analysis_name#'.' - '.'#species_name#',
     },
  };
  
  
  ## SINGLECELL: convert bam file to cram
  push @pipeline, {
    -logic_name  => 'convert_cellranger_bam_to_cram',
    -module      => 'ehive.runnable.process.alignment.ConvertBamToCram',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
        'bam_file'        => '#bam_file#',
        'samtools_exe'    => $self->o('samtools_exe'),
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
  
  
  ## SINGLECELL: upload cram file to irods server
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
      'analysis_dir'  => $self->o('analysis_dir'),
      'dir_path_list' => ['#analysis_dir#','#sample_igf_id#','#experiment_igf_id#','#analysis_name#'],
      'file_tag'      => '#sample_igf_id#'.' - '.'#experiment_igf_id#'.' - '.'#analysis_name#'.' - '.'#species_name#',
     },
  };
  
  
  ## SINGLECELL: scanpy report
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
  
  
  ## SINGLECELL: copy scanpy report to remote
  push @pipeline, {
      -logic_name   => 'copy_scanpy_report_to_remote',
      -module       => 'ehive.runnable.process.alignment.CopyAnalysisFilesToRemote',
      -language     => 'python3',
      -meadow_type  => 'PBSPro',
      -rc_name      => '1Gb',
      -analysis_capacity => 2,
      -parameters  => {
        'analysis_dir'        => $self->o('analysis_dir'),
        'dir_labels'          => ['#analysis_dir#','#sample_igf_id#'],
        'file_list'           => ['#output_report#'],
        'remote_user'         => $self->o('seqrun_user'),
        'remote_host'         => $self->o('remote_host'),
        'remote_project_path' => $self->o('remote_project_path'),
        },
  };
  
  
  ## SINGLECELL: picard alignment summary metrics
  push @pipeline, {
    -logic_name  => 'picard_aln_summary_for_cellranger',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => ['#bam_file#'],
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'CollectAlignmentSummaryMetrics',
      'base_work_dir'  => $self->o('base_work_dir'),
      'reference_type' => $self->o('reference_fasta_type'),
     },
    -flow_into   => {
        1 => ['picard_base_dist_summary_for_cellranger'],
      },
  };
  
  
  ## SINGLECELL: picard base distribution summary metrics
  push @pipeline, {
    -logic_name  => 'picard_base_dist_summary_for_cellranger',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => ['#bam_file#'],
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'CollectBaseDistributionByCycle',
      'base_work_dir'  => $self->o('base_work_dir'),
      'reference_type' => $self->o('reference_fasta_type'),
     },
    -flow_into   => {
        1 => ['picard_gc_bias_summary_for_cellranger'],
      },
  };
  
  
  ## SINGLECELL: picard gc bias summary metrics
  push @pipeline, {
    -logic_name  => 'picard_gc_bias_summary_for_cellranger',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => ['#bam_file#'],
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'CollectGcBiasMetrics',
      'base_work_dir'  => $self->o('base_work_dir'),
      'reference_type' => $self->o('reference_fasta_type'),
     },
    -flow_into   => {
        1 => ['picard_qual_dist_summary_for_cellranger'],
      },
  };
  
  
  ## SINGLECELL: picard quality distribution summary metrics
  push @pipeline, {
    -logic_name  => 'picard_qual_dist_summary_for_cellranger',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => ['#bam_file#'],
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'QualityScoreDistribution',
      'base_work_dir'  => $self->o('base_work_dir'),
      'reference_type' => $self->o('reference_fasta_type'),
     },
    -flow_into   => {
        1 => ['picard_rna_metrics_summary_for_cellranger'],
      },
  };
  
  
  ## SINGLECELL: picard rna metrics summary metrics
  push @pipeline, {
    -logic_name  => 'picard_rna_metrics_summary_for_cellranger',
    -module      => 'ehive.runnable.process.alignment.RunPicard',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '4Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'    => ['#bam_file#'],
      'java_exe'       => $self->o('java_exe'),
      'java_param'     => $self->o('java_param'),
      'picard_jar'     => $self->o('picard_jar'),
      'picard_command' => 'CollectRnaSeqMetrics',
      'base_work_dir'  => $self->o('base_work_dir'),
      'reference_type' => $self->o('reference_fasta_type'),
      'reference_refFlat' => $self->o('reference_refFlat'),
     },
    -flow_into   => {
        1 => ['samtools_flagstat_summary_for_cellranger'],
      },
  };
  
  
  ## SINGLECELL: samtools stats metrics
  push @pipeline, {
    -logic_name  => 'samtools_flagstat_summary_for_cellranger',
    -module      => 'ehive.runnable.process.alignment.RunSamtools',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'      => ['#bam_file#'],
      'samtools_exe'     => $self->o('samtools_exe'),
      'samtools_command' => 'stats',
      'base_work_dir'    => $self->o('base_work_dir'),
      'reference_type'   => $self->o('reference_fasta_type'),
      'threads'          => $self->o('samtools_threads'),
     },
    -flow_into   => {
        1 => ['samtools_idxstat_summary_for_cellranger'],
      },
  };
  
  
  ## SINGLECELL: samtools idxstat metrics
  push @pipeline, {
    -logic_name  => 'samtools_idxstat_summary_for_cellranger',
    -module      => 'ehive.runnable.process.alignment.RunSamtools',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'      => ['#bam_file#'],
      'samtools_exe'     => $self->o('samtools_exe'),
      'samtools_command' => 'idxstats',
      'base_work_dir'    => $self->o('base_work_dir'),
      'reference_type'   => $self->o('reference_fasta_type'),
     },
    -flow_into   => {
        1 => ['multiqc_report_for_cellranger'],
      },
  };
  
  
  ## SINGLECELL: multiqc report building
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
      'analysis_name'    => $self->o('multiqc_analysis'),
      'tag_name'         => '#species_name#',
      'multiqc_exe'      => $self->o('multiqc_exe'),
      'multiqc_options'  => $self->o('multiqc_options'),
     },
    -flow_into   => {
        1 => ['copy_sample_multiqc_for_singlecell_to_remote'],
      },
  };
  
  
  ## SINGLECELL: copy multiqc to remote
  push @pipeline, {
      -logic_name   => 'copy_sample_multiqc_for_singlecell_to_remote',
      -module       => 'ehive.runnable.process.alignment.CopyAnalysisFilesToRemote',
      -language     => 'python3',
      -meadow_type  => 'PBSPro',
      -rc_name      => '1Gb',
      -analysis_capacity => 2,
      -parameters  => {
        'analysis_dir'        => $self->o('analysis_dir'),
        'dir_labels'          => ['#analysis_dir#','#sample_igf_id#'],
        'file_list'           => ['#multiqc_html#'],
        'remote_user'         => $self->o('seqrun_user'),
        'remote_host'         => $self->o('remote_host'),
        'remote_project_path' => $self->o('remote_project_path'),
        },
      -flow_into    => {
          1 => ['mark_experiment_finished'],
      },
  };
  
  
  #############################  SINGLECELL END   ##############################
  
  
  ## GENERIC: mark experiment as done
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
        'rna_source'    => $self->o('rna_source'),
        
        },
       -flow_into    => {
          1 => ['update_project_analysis'],
      },
  };
  
  
  return \@pipeline;
}

1;

