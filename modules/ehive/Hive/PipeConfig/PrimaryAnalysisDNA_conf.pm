=head1 NAME
    ehive::Hive::PipeConfig::PrimaryAnalysisDNA_conf
=cut


package ehive::Hive::PipeConfig::PrimaryAnalysisCombinedDNA_conf;

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
    'patterned_flow_cell_list'     => ['NEXTSEQ','HISEQ4000'],
    'rna_source'          => 'TRANSCRIPTOMIC',
    ## Irods
    'irods_exe_dir'       => undef,
    ## Java
    'java_exe'            => undef,
    'java_param'          => '-Xmx4g',
    ## Picard
    'picard_jar'                   => undef,
    'illumina_platform_name'       => 'ILLUMINA',
    ## MultiQC
    'multiqc_analysis'    => 'multiqc',
    'multiqc_exe'         => undef,
    'multiqc_options'     => '{"--zip-data-dir" : ""}',
    'multiqc_type'        => 'MULTIQC_HTML',
    ## Ref genome
    'reference_fasta_type'=> 'GENOME_FASTA',
    'reference_refFlat'   => 'GENE_REFFLAT',
    'reference_gtf_type'  => 'GENE_GTF',
    'two_bit_genome_type' => 'GENOME_TWOBIT_URI',
    ## Fetch fastq
    'fastq_collection_type'        => undef,
    'fastq_collection_table'       => undef,
    ## Fastp adapter trimming
    'fastp_exe'            => undef,
    'fastp_options_list'   => ['--qualified_quality_phred=15','--length_required=15'],
    'fastp_run_thread'     => 4,
    'split_by_lines_count' => 5000000,
    'fastp_analysis_name'  => 'fastp',
    'fastp_html_collection_type'   => 'FASTP_REPORT',
    'fastp_collection_table'       => 'run',
    ## Samtools
    'samtools_exe'         => undef,
    'samtools_threads'     => 4,
    ## STAR alignment
    'star_exe'             => undef,
    'star_reference_type'  => 'TRANSCRIPTOME_STAR',
    'star_patameters'      => '{"--outFilterMultimapNmax":"20","--alignSJoverhangMin":"8","--alignSJDBoverhangMin":"1","--outFilterMismatchNmax":"999","--outFilterMismatchNoverReadLmax":"0.04","--alignIntronMin":"20","--alignIntronMax":"1000000","--alignMatesGapMax":"1000000","--outSAMattributes":"NH HI AS NM MD","--limitBAMsortRAM":"12000000000"}',
    'star_run_thread'      => 8,
    'star_two_pass_mode'   => 1,
    'star_analysis_name'   => undef,
    'star_multiqc_type'    => undef,
    'bedGraphToBigWig_path'        => undef,
    'star_collection_table'        => undef,
    'star_genomic_cram_type'       => undef,
    'star_bw_collection_type'      => undef,
    ## BWA
    'bwa_exe'              => undef,
    'bwa_reference_type'   => undef,
    'bwa_run_thread'       => undef,
    'bwa_parameters'       => '{-"M":""}',
    ## RSEM
    'rsem_exe_dir'         => undef,
    'rsem_reference_type'  => 'TRANSCRIPTOME_RSEM',
    'rsem_analysis_name'   => 'rsem',
    'rsem_threads'         => 8,
    'rsem_memory_limit'    => 4000,
    'rsem_analysis_name'   => undef,
    'rsem_collection_type' => undef,
    'rsem_collection_table'        => undef,
    ## Cellranger
    'cellranger_exe'       => undef,
    'cellranger_param'     => '{"--nopreflight":"","--disable-ui":"","--jobmode":"pbspro","--localcores":"1","--localmem":"1","--mempercore":"4","--maxjobs":"20"}',
    'cellranger_timeout'   => 43200,
    'cellranger_collection_table'  => 'experiment',
    'cellranger_analysis_name'     => 'cellranger_count',
    ## Scanpy
    'scanpy_type'          => 'SCANPY_RESULTS',
    'scanpy_report_template'       => undef,
    ## Demultiplexing pipeline
    'demultiplexing_pipeline_name' => undef,
    ## Remote dir settings
    'seqrun_user'          => undef,
    'remote_host'          => undef,
    'remote_project_path'  => undef,
    'analysis_dir'         => 'analysis',
    ## Genome browser
    'genome_browser_template_file' => undef,
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
      'rna_source'    => $self->o('rna_source'),
      'dna_source'    => $self->o('genomic_source'),
    },
    -flow_into => {
        2 => WHEN('#library_source# eq #dna_source#' => ['run_factory_for_dnaseq'],
                       ELSE ['mark_experiment_finished']),
    },
  };
  
  
  ## run factory for genomic data
  push @pipeline, {
    -logic_name  => 'run_factory_for_dnaseq',
    -module      => 'ehive.runnable.jobfactory.alignment.RunFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -analysis_capacity => 2,
    -flow_into   => {
        '2->A' => ['fetch_fastq_for_run'],
        'A->1' => ['process_bwa_bams'],
      },
  };
  
  
  ## fetch fastq files for a run
  push @pipeline, {
    -logic_name  => 'fetch_fastq_for_run',
    -module      => 'ehive.runnable.process.alignment.FetchFastqForRun',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -rc_name     => '1Gb',
    -analysis_capacity => 5,
    -parameters  => {
      'fastq_collection_type'  => $self->o('fastq_collection_type'),
      'fastq_collection_table' => $self->o('fastq_collection_table'),
    },
    -flow_into   => {
        1 => ['adapter_trim_without_fastq_split'],
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
      'fastp_options_list'   => $self->o('fastp_options_list'),
      'split_by_lines_count' => $self->o('split_by_lines_count'),
      'run_thread'           => $self->o('fastp_run_thread'),
      'base_work_dir'        => $self->o('base_work_dir'),
      'fastp_exe'            => $self->o('fastp_exe'),
      'input_fastq_list'     => '#fastq_files_list#',
    },
    -flow_into   => {
        '1' => ['load_fastp_report'],
      },
  };
  
  
  ## collect fastp report
  push @pipeline, {
    -logic_name  => 'load_fastp_report',
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
     },
    -flow_into   => {
        1 => ['run_bwa'],
      },
  };
  
  
  ## run bwa alignment
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
          1 => ['picard_add_rg_tag_to_genomic_bam'],
    },
  };
  
  
  ## picard add rg tag to run star genomic bam
  push @pipeline, {
    -logic_name  => 'picard_add_rg_tag_to_genomic_bam',
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
         'RGCN'       => 'Imperial Genomics Facility',
         'SORT_ORDER' => 'coordinate',
         },
     },
    -flow_into => {
          1 => [ '?accu_name=bwa_aligned_genomic_bam&accu_address={experiment_igf_id}{seed_date_stamp}[]&accu_input_variable=analysis_files' ],
     },
  };
  
  
  ## process bwa bams
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
  
  
  ## collect star genomic bam
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
        1 => {'picard_merge_and_mark_dup_bwa_genomic_bam' => {'bwa_genomic_bams' => '#exp_chunk_list#'}},
      },
  };
  
  
  ## picard merge and mark duplicate genomic bam
  push @pipeline, {
    -logic_name  => 'picard_merge_and_mark_dup_bwa_genomic_bam',
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
      'SORT_ORDER'     => 'coordinate',
     },
    -flow_into => {
          1 => { 'convert_bwa_genomic_bam_to_cram' => {'merged_bwa_genomic_bams' => '#bam_files#',
                                                       'analysis_files' => '#analysis_files#'}},
     },
  };
  
  
  ## convert genomic bam to cram
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
  
  
  ## copy star genomic cram to irods
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
  
  
  ## picard alignment summary metrics for bwa
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
     },
    -flow_into   => {
        1 => ['picard_base_dist_summary_for_bwa'],
      },
  };
  
  
  ## picard base distribution summary metrics
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
     },
    -flow_into   => {
        1 => ['picard_gc_bias_summary_for_bwa'],
      },
  };
  
  
  ## picard base distribution summary metrics
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
     },
    -flow_into   => {
        1 => ['picard_gc_bias_summary_for_bwa'],
      },
  };
  
  
  ## picard gc bias summary metrics
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
     },
    -flow_into   => {
        1 => ['picard_qual_dist_summary_for_bwa'],
      },
  };
  
  
  ## picard quality distribution summary metrics
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
     },
    -flow_into   => {
        1 => ['samtools_flagstat_summary_for_bwa'],
      },
  };
  
  
  ## samtools flagstat metrics
  push @pipeline, {
    -logic_name  => 'samtools_flagstat_summary_for_bwa',
    -module      => 'ehive.runnable.process.alignment.RunSamtools',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb4t',
    -analysis_capacity => 2,
    -parameters  => {
      'input_files'      => '#merged_bwa_genomic_bams#',
      'samtools_command' => 'flagstat',
      'base_work_dir'    => $self->o('base_work_dir'),
      'reference_type'   => $self->o('reference_fasta_type'),
      'samtools_exe'     => $self->o('samtools_exe'),
      'threads'          => $self->o('samtools_threads'),
     },
    -flow_into   => {
        1 => ['samtools_idxstat_summary_for_bwa'],
      },
  };
  
  
  ## samtools idxstat metrics
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
      'base_work_dir'    => $self->o('base_work_dir'),
      'samtools_exe'     => $self->o('samtools_exe'),
      'reference_type'   => $self->o('reference_fasta_type'),
     },
    -flow_into   => {
        1 => ['multiqc_report_for_bwa'],
      },
  };
  
  
  ## multiqc report building
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
      'collection_type'  => $self->o('bwa_multiqc_type'),
      'collection_table' => $self->o('bwa_collection_table'),
      'analysis_name'    => $self->o('multiqc_analysis'),
      'tag_name'         => '#species_name#',
      'multiqc_exe'      => $self->o('multiqc_exe'),
      'multiqc_options'  => $self->o('multiqc_options'),
     },
    -flow_into   => {
        1 => ['copy_bwa_multiqc_to_remote'],
      },
  };
  
  
  ## copy multiqc to remote
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
        'rna_source'    => $self->o('rna_source'),
        
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
      -flow_into    => {
          1 => ['update_project_status'],
      },
  };
  
  
  ## update status page
  push @pipeline, {
      -logic_name   => 'update_project_status',
      -module       => 'ehive.runnable.process.UpdateProjectStatus',
      -language     => 'python3',
      -meadow_type  => 'PBSPro',
      -rc_name      => '1Gb',
      -analysis_capacity => 1,
      -parameters   => {
        'remote_project_path'  => $self->o('remote_project_path'),
        'remote_user'          => $self->o('seqrun_user'),
        'remote_host'          => $self->o('remote_host'),
        'demultiplexing_pipeline_name' => $self->o('demultiplexing_pipeline_name'),
        'analysis_pipeline_name'       => $self->o('pipeline_name'),
      },
  };
  
  return \@pipeline;
}

1;