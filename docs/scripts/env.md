# Pipeline source code and environment

## Checkout codebases from Github

* Pipeline config
  <p>```
  git clone https://github.com/imperial-genomics-facility/igf-pipelines.git</p>
  ```</p>

* Python library for IGF pipelines
  <p>```
  git clone https://github.com/imperial-genomics-facility/data-management-python.git
  ```</p>

* Ensembl ehive pipeline
  <p>```
  git clone https://github.com/Ensembl/ensembl-hive.git
  ```</p>

* PBSpro meadow interface for Ensembl ehive
  <p>```
  git clone https://github.com/Ensembl/ensembl-hive-pbspro.git
  ```</p>

## Install required perl packages

Use system perl or install specific perl version using [perlbrew](https://perlbrew.pl/). Check for required perl modules after installation using the following command
  <p>```
  perl -e 'use DBI;use DBD::mysql;use DBD::SQLite;IPC::Run;'
  ```</p>

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
  bash Miniconda3-latest-Linux-x86_64.sh
</code></pre>

* Create an conda environment using the configuration file
  <p>```
  conda env create --name pipeline-env --file igf-pipelines/environment.yaml
  ```</p>

## Test conda environment

Activate the conda environment using the following command
  <p>```
  source activate pipeline-env
  ```</p>

Deactivate conda env
  <p>```
  source deactivate
  ```</p>

## Create environment configuration file for pipeline

Create a new file, e.g. `env.sh` and 

 *  Add conda activation command
  <p>```
  source activate pipeline-env
  ```</p>
  
 * Add PERL5LIB variables
  <pre><code>  export PERL5LIB=/path/ensembl-hive/modules:/path/ensembl-hive-pbspro/modules:/path/igf-pipelines/modules:${PERL5LIB}
  </code></pre>
  
  * Add PATH variable
  <pre><code>  export PATH=/path/ensembl-hive/scripts:${PATH}
  </code></pre>
  
  * Add PYTHONPATH variables
  <pre><code>  export PYTHONPATH=/path/data-management-python/:/path/ensembl-hive/wrappers/python3:${PYTHONPATH}
  </code></pre>

Then activate environment for pipeline using `source env.sh` command.
