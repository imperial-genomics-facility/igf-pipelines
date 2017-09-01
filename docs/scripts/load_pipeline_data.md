# Load pipeline data to db

A script for loading pipeline data to database

## Usage
<pre><code>load_pipeline_data.py [-h] -p PIPELINE_DATA [-u] -d DBCONFIG_PATH -s SLACK_CONFIG</pre></code>

## Parameters
<pre><code>  -h /--help     show this help message and exit
  -p /--pipeline_data    Pipeline data json file
  -u /--update           Update existing platform data, default: False
  -d /--dbconfig_path    Database configuration json file
  -s /--slack_config     Slack configuration json file
</code></pre>

## Pipeline data json file
<pre><code>{ 
 "pipeline_name":"demultiplexing_fastq",
 "pipeline_db":"sqlite:////data/bcl2fastq.db", 
 "pipeline_init_conf":{"input_dir":"data/seqrun_dir/", "output_dir":"data"}, 
 "pipeline_run_conf":{"output_dir":"/path/output_dir2"},
}</code></pre>
