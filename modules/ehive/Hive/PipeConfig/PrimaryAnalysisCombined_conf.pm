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
  
  return \@pipeline;
}

1;

