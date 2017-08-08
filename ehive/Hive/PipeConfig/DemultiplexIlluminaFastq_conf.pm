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
    -logic_name => 'find_new_sequencing_runs',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
    -flow_into => {
      '2->A' => [ 'find_customer_information' ],
      'A->1' => [ 'mark_sequencing_complete' ],
    },
  };
  
  push @pipeline, {
    -logic_name => 'find_customer_information',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => [ 'setup_customer_account' ],
    },
  };
  
  push @pipeline, {
    -logic_name => 'setup_customer_account',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => [ 'file_md5_pre_transfer' ],
    },
  };
  
  push @pipeline, {
    -logic_name => 'file_md5_pre_transfer',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => [ 'find_seqrun_files' ],
    },
  };
  
  push @pipeline, {
    -logic_name => 'find_seqrun_files',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
    -flow_into => {
      '2->A' => ['transfer_seqrun_file'],
      'A->1' => ['find_seqrun_lanes'],
    },
  };
  
  push @pipeline, {
    -logic_name => 'transfer_seqrun_file',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => [ 'file_md5_post_transfer' ],
    },
  };
    
  push @pipeline, {
    -logic_name => 'file_md5_post_transfer',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => [ '?accu_name=file_md5&accu_address={file}&accu_input_variable=md5_value' ],
    }, 
  };
  
  push @pipeline, {
    -logic_name => 'find_seqrun_lanes',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
    -flow_into => {
      '2->A' => [ 'prepare_samplesheet' ],
      'A->1' => [ 'generate_fastqc_report' ],
    },
  };
  
  push @pipeline, {
    -logic_name => 'prepare_samplesheet',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => { 'run_bcl2fastq' => {'samplesheet'=>'#samplesheet#', 'lane' => '#lane#'}},
    },
  };
  
  push @pipeline, {
    -logic_name => 'run_bcl2fastq',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => [ 'check_undetermined_barcodes' ],
    },
  };
  
  push @pipeline, {
    -logic_name => 'check_undetermined_barcodes',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => {'calculate_exp_and_run_from_samplesheet' => {'fastq_file' => '#fastq_file#' }},
    },
  };
  
  push @pipeline, {
    -logic_name => 'calculate_exp_and_run_from_samplesheet',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => [ 'collect_fastq' ],
    },
  };
  
  push @pipeline, {
    -logic_name => 'collect_fastq',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => {'upload_demultiplexing_report' => {'fastq_files' => '#fastq_files#', 'lane' => '#lane#', 'report' => '#report#', 'undetermined_fastqs' => '#undetermined_fastqs#' }},
    },
  };
  
  push @pipeline, {
    -logic_name => 'upload_demultiplexing_report',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => {'upload_to_irods' => {'fastq_files' => '#fastq_files#', 'lane' => '#lane#' }},
    },
  };
  
  push @pipeline, {
    -logic_name => 'upload_to_irods',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      '2->A' => {'find_fastq_for_lane' => {'fastq_files' => '#fastq_files#',  }},
      'A->1' => {'find_undetermined_fastq_for_lane' => {'undetermined_fastqs' => '#undetermined_fastqs#'  }},
    },
  };
  
  push @pipeline, {
    -logic_name => 'find_fastq_for_lane',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
    -flow_into => {
      '2->A' => [ 'run_fastqc' ],
      'A->1' => [ 'run_multiqc' ],
    },
  };
    
  push @pipeline, {
    -logic_name => 'run_fastqc',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => {'upload_fastqc' => {'fastqc_file' => '#fastqc_file#',  }},
    },
  };
  
  push @pipeline, {
    -logic_name => 'upload_fastqc',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => [ '?accu_name=fastq_file_list&accu_address={fastq_file}&accu_input_variable=fastqc_file' ],
    },
  };
  
  push @pipeline, {
    -logic_name => 'run_multiqc',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => {'upload_multiqc' => {'multi_qc_report' => '#multi_qc_report#',  }},
    },
  };
  
  push @pipeline, {
    -logic_name => 'upload_multiqc',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
  };
  
  push @pipeline, {
    -logic_name => 'generate_fastqc_report',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => [ 'upload_fastqc_report' ],
    },
  };
  
  push @pipeline, {
    -logic_name => 'upload_fastqc_report',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
  };
  
  push @pipeline, {
    -logic_name => 'find_undetermined_fastq_for_lane',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
    -flow_into => {
      '2->A' => [ 'run_undetermined_fastqc' ],
      'A->1' => [ 'add_undetermined_barcodes' ],
    },
  };
  
  push @pipeline, {
    -logic_name => 'run_undetermined_fastqc',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -flow_into => {
      1 => {'upload_undetermined_fastqc' => {'undetermined_fastqc_file' => '#undetermined_fastqc_file#',  }},
    },
  };
  
  push @pipeline, {
    -logic_name => 'upload_undetermined_fastqc',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
  };
  
  push @pipeline, {
    -logic_name => 'add_undetermined_barcodes',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
  };
  
  push @pipeline, {
    -logic_name => 'mark_sequencing_complete',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
  };
  
 
  return \@pipeline
}

1;
