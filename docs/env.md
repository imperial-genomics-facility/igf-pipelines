# Pipeline source code and environment

## Checkout codebases from Github

* Pipeline config
<pre><code>
  git clone https://github.com/imperial-genomics-facility/igf-pipelines.git

</code></pre>

* Python library for IGF pipelines
<pre><code>
  git clone https://github.com/imperial-genomics-facility/data-management-python.git

</code></pre>

* Ensembl ehive pipeline
<pre><code>
  git clone https://github.com/Ensembl/ensembl-hive.git

</code></pre>

* PBSpro meadow interface for Ensembl ehive
<pre><code>
  git clone https://github.com/Ensembl/ensembl-hive-pbspro.git

</code></pre>

## Install required perl and python packages
Setup environment using conda

* Download and setup miniconda
<pre><code>
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
  
  bash Miniconda3-latest-Linux-x86_64.sh

</code></pre>

* Create an conda environment using the configuration file
<pre><code>
  conda env create --name pipeline-env --file igf-pipelines/environment.yaml

</code></pre>

## Set additional PERL5LIB and PYTHONPATH variables
Add the following lines to the `~/.bashrc` file
<pre><code>
  export PERL5LIB=/path/ensembl-hive/modules:/path/ensembl-hive-pbspro/modules:/path/igf-pipelines/modules:${PERL5LIB}
  export PATH=/path/ensembl-hive/scripts:${PATH}
  export PYTHONPATH=/path/data-management-python/:/path/ensembl-hive/wrappers/python3:${PYTHONPATH}

</code></pre>

## Activate conda environment
Activate the conda environment using the following command
<pre><code>
  source activate pipeline-env

</code></pre>

## Install missing perl modules
Check for required perl modules after activating the environment using the following command
<pre><code>
  perl -e 'use DBI;use DBD::mysql;use DBD::SQLite;IPC::Run;'

</code></pre>

Install the missing packages using `cpanm` if they are not already present.
<pre><code>
  cpanm DBI
  cpanm DBD::mysql
  cpanm DBD::SQLite
  cpanm IPC::Run

</code></pre>
