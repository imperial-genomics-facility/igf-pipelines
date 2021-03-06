=head1 NAME
    ehive::Hive::PipeConfig::IGFBasePipe_conf

=cut

package ehive::Hive::PipeConfig::IGFBasePipe_conf;


use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;


sub default_options {
  my $self = shift;
  
  return {
      %{$self->SUPER::default_options},
      pbs_queue        => 'pqcgi',
      user_email       => '/dev/null',
      log_slack        => 1,
      log_asana        => 1,
      dbconfig         => undef,
      asana_config     => undef,
      slack_config     => undef,
      asana_project_id => undef,
   };
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},            # here we inherit anything from the base class
        'dbconfig'     => $self->o('dbconfig'),
        'log_slack'    => $self->o('log_slack'),
        'slack_config' => $self->o('slack_config'),
        'log_asana'    => $self->o('log_asana'),
        'asana_config' => $self->o('asana_config'),
        'asana_project_id' => $self->o('asana_project_id'),
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
            %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            '500Mb2t'  => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=2:mem=500mb' },
            '500Mb'    => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=1:mem=500mb' },
            '1Gb'      => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=1:mem=1gb' },
            '1Gb2t'    => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=2:mem=1gb' },
            '1Gb4t'    => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=4:mem=1gb' },
            '2Gb'      => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=1:mem=2gb' },
            '2Gb72hr'  => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=72:00:00 -l select=1:ncpus=1:mem=2gb' },
            '2Gb2t'    => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=2:mem=2gb' },
            '2Gb4t'    => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=4:mem=2gb' },
            '2Gb8t'    => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=8:mem=2gb' },
            '4Gb'      => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=1:mem=4gb' },
            '4Gb2t'    => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=2:mem=4gb' },
            '4Gb4t'    => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=4:mem=4gb' },
            '4Gb8t'    => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=8:mem=4gb' },
            '4Gb16t'   => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=16:mem=4gb' },
            '8Gb'      => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=1:mem=8gb' },
            '8Gb2t'    => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=2:mem=8gb' },
            '8Gb4t'    => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=4:mem=8gb' },
            '8Gb8t'    => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=8:mem=8gb' },
            '8Gb16t'   => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=16:mem=8gb' },
            '12Gb'     => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=1:mem=12gb' },
            '12Gb4t'   => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=4:mem=12gb' },
            '12Gb8t'   => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=8:mem=12gb' },
            '16Gb'     => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=1:mem=16gb' },
            '16Gb4t'   => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=4:mem=16gb' },
            '16Gb8t'   => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=8:mem=16gb' },
            '16Gb16t'  => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=16:mem=16gb' },
            '32Gb'     => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=1:mem=32gb' },
            '32Gb4t'   => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=4:mem=32gb' },
            '32Gb8t'   => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=8:mem=32gb' },
            '32Gb16t'  => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=16:mem=32gb' },
            '42Gb4t'   => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=4:mem=42gb' },
            '42Gb8t'   => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=8:mem=42gb', },
            '42Gb16t'  => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=24:00:00 -l select=1:ncpus=16:mem=42gb', },
            '64Gb8t72hr'   => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=72:00:00 -l select=1:ncpus=8:mem=64gb', },
            '64Gb16t72hr'  => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=72:00:00 -l select=1:ncpus=16:mem=64gb', },
            '128Gb16t72hr' => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -k n -m n -o /dev/null -e /dev/null -l walltime=72:00:00 -l select=1:ncpus=16:mem=128gb', },
    };
}


sub hive_meta_table {
  my ($self) = @_;
  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack'  => 1,           # switching on parameter propagation by default
    'hive_auto_rebalance_semaphores' => 1,
  };
}

1;
