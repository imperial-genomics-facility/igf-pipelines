=head1 NAME
    ehive::Hive::PipeConfig::DemultiplexIlluminaFastq_conf

=cut

package ehive::Hive::PipeConfig::DemultiplexIlluminaFastq_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('ehive::Hive::PipeConfig::IGFBasePipe_conf');


sub pipeline_analyses {
  my ($self) = @_;
  my @pipeline;
  
  push @pipeline, {
    -logic_name  => 'find_new_sequencing_runs',
    -module      => 'ehive.runnable.jobfactory.PipeseedFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -flow_into => {
      2 => ['seqrun_transfer_factory'],    
    },
  };
  push @pipeline, {
    -logic_name  => 'seqrun_transfer_factory',
    -module      => 'ehive.runnable.jobfactory.SeqrunFileFactory',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
    -flow_into => {
      2 => ['transfer_seqrun_file'],    
    },
  };
  push @pipeline, {
    -logic_name  => 'transfer_seqrun_file',
    -module      => 'ehive.runnable.proccess.TransferAndCheckRemoteFile',
    -language    => 'python3',
    -meadow_type => 'LOCAL',
  };
  
 
  return \@pipeline
}

1;
