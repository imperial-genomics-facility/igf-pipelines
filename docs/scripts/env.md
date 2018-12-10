# Pipeline source code and environment

## Checkout codebases from Github

* Pipeline config
  <pre><code>  git clone https://github.com/imperial-genomics-facility/igf-pipelines.git  </code></pre>

* Python library for IGF pipelines
  <pre><code>  git clone https://github.com/imperial-genomics-facility/data-management-python.git  </code></pre>

* Ensembl ehive pipeline
  <pre><code>  git clone https://github.com/Ensembl/ensembl-hive.git  </code></pre>

* PBSpro meadow interface for Ensembl ehive
  <pre><code>  git clone https://github.com/Ensembl/ensembl-hive-pbspro.git  </code></pre>

## Install required perl packages

Use system perl or install specific perl version using [perlbrew](https://perlbrew.pl/). Check for required perl modules after installation using the following command
  <pre><code>  perl -e 'use DBI;use DBD::mysql;use DBD::SQLite;IPC::Run;'  </code></pre>

Install the missing packages using `cpan` (for system perl) or `cpanm` (for perlbrew) if they are not already present.
<pre><code>  cpan install DBI
  cpan install DBD::mysql
  cpan install DBD::SQLite
  cpan install IPC::Run  
</code></pre>


## Install required python packages

Setup environment using conda

* Download and setup miniconda
  <pre><code>  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh  </code></pre>

* Create an conda environment using the configuration file
  <pre><code>  conda env create --name pipeline-env --file igf-pipelines/environment.yaml  </code></pre>

## Test conda environment

Activate the conda environment using the following command
  <pre><code>  source activate pipeline-env  </code></pre>

Deactivate conda env
  <pre><code>  source deactivate  </code></pre>

## Create environment configuration file for pipeline

Create a new file, e.g. `env.sh` and 

 *  Add conda activation command
  <pre><code>  source activate pipeline-env  </code></pre>
  
 * Add PERL5LIB variables
  <pre><code>  export PERL5LIB=/path/ensembl-hive/modules:/path/ensembl-hive-pbspro/modules:/path/igf-pipelines/modules:${PERL5LIB}  </code></pre>
  
  * Add PATH variable
  <pre><code>  export PATH=/path/ensembl-hive/scripts:${PATH}  </code></pre>
  
  * Add PYTHONPATH variables
  <pre><code>  export PYTHONPATH=/path/data-management-python/:/path/ensembl-hive/wrappers/python3:${PYTHONPATH}  </code></pre>

Then activate environment for pipeline using `source env.sh` command.
