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

## Install required perl and python packages
Setup environment using conda

* Download and setup miniconda
 <pre><code>  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
</code></pre>

* Create an conda environment using the configuration file
  <p>```
  conda env create --name pipeline-env --file igf-pipelines/environment.yaml
  ```</p>


## Set additional PERL5LIB and PYTHONPATH variables
Add the following lines to the `~/.bashrc` file
<pre><code>  export PERL5LIB=/path/ensembl-hive/modules:/path/ensembl-hive-pbspro/modules:/path/igf-pipelines/modules:${PERL5LIB}
  export PATH=/path/ensembl-hive/scripts:${PATH}
  export PYTHONPATH=/path/data-management-python/:/path/ensembl-hive/wrappers/python3:${PYTHONPATH}
</code></pre>

## Activate conda environment
Activate the conda environment using the following command
  <p>```
  source activate pipeline-env
  ```</p>

## Install missing perl modules
Check for required perl modules after activating the environment using the following command
  ```
  perl -e 'use DBI;use DBD::mysql;use DBD::SQLite;IPC::Run;'
  ```

Install the missing packages using `cpanm` if they are not already present.
<pre><code>  cpanm DBI
  cpanm DBD::mysql
  cpanm DBD::SQLite
  cpanm IPC::Run
</code></pre>
