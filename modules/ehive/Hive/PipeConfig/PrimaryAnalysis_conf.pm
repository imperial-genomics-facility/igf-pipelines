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
    'reference_fasta_type'        => 'GENOME_FASTA',
    'cellranger_collection_table' => 'experiment',
    'cellranger_analysis_name'    => 'cellranger_count',
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
    -logic_name  => 'find_new_sequencing_runs',
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
      },
    -flow_into   => {
        1 => ['convert_bam_to_cram'],  
      },
  };
  
  push @pipeline, {
    -logic_name  => 'convert_bam_to_cram',
    -module      => 'ehive.runnable.process.alignment.ConvertBamToCram',
    -language    => 'python3',
    -meadow_type => 'PBSPro',
    -rc_name     => '2Gb',
    -analysis_capacity => 1,
    -parameters  => {
        'bam_file'        => '#bam_file#',
        'base_result_dir' => $self->o('base_results_dir'),
        'collection_name' => '#experiment_igf_id#',
        'collection_type' => $self->o('cram_type'),
        'collection_table'=> $self->o('cellranger_collection_table'),
        'analysis_name'   => $self->o('cellranger_analysis_name'),
        'tag_name'        => '#species_name#',
        'reference_type'  => $self->o('reference_fasta_type'),
     },
    -flow_into   => {
        1 => ['mark_experiment_finished'],  
      },
  };
  
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