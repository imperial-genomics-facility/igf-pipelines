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
  
  return \@pipeline;
}

1;