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
      pbs_queue => 'pqcgi',
      user_email => '/dev/null',
   };
}


sub resource_classes {
    my ($self) = @_;
    return {
            %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            '500Mb' => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -M /dev/null -l walltime=72:00:00 -l select=1:ncpus=1:mem=500mb:tmpspace=10gb' },
            '500Mb4t' => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -M /dev/null -l walltime=72:00:00 -l select=1:ncpus=4:mem=500mb:tmpspace=10gb' },
            '1Gb' => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -M /dev/null -l walltime=72:00:00 -l select=1:ncpus=1:mem=1gb:tmpspace=10gb' },
            '1Gb4t' => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -M /dev/null -l walltime=72:00:00 -l select=1:ncpus=4:mem=1gb:tmpspace=10gb' },
            '2Gb' => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -M /dev/null -l walltime=72:00:00 -l select=1:ncpus=1:mem=2gb:tmpspace=10gb' },
            '2Gb4t' => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -M /dev/null -l walltime=72:00:00 -l select=1:ncpus=4:mem=2gb:tmpspace=10gb' },
            '2GbDebug' => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -M '.$self->o('user_email').' -l walltime=72:00:00 -l select=1:ncpus=1:mem=2gb:tmpspace=10gb' },
            '4Gb' => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -M /dev/null -l walltime=72:00:00 -l select=1:ncpus=1:mem=4gb:tmpspace=10gb' },
            '4Gb4t' => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -M /dev/null -l walltime=72:00:00 -l select=1:ncpus=4:mem=4gb:tmpspace=10gb' },
            '4GbDebug' => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -M '.$self->o('user_email').' -l walltime=72:00:00 -l select=1:ncpus=1:mem=4gb:tmpspace=10gb' },
            '8Gb' => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -M /dev/null -l walltime=72:00:00 -l select=1:ncpus=1:mem=8gb:tmpspace=10gb' },
            '8Gb4t' => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -M /dev/null -l walltime=72:00:00 -l select=1:ncpus=4:mem=8gb:tmpspace=10gb' },
            '8Gb16t' => { 'PBSPro' => '-q '.$self->o('pbs_queue').' -M /dev/null -l walltime=72:00:00 -l select=1:ncpus=16:mem=8gb:tmpspace=10gb' },
            
    };
}


sub hive_meta_table {
  my ($self) = @_;
  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack'  => 1,           # switching on parameter propagation by default
  };
}

1;
