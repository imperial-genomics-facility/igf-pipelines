# Find new sequencing run in server

## Description
This script is for identifiying new finished sequencing runs in the server. It will prepare a file list including the md5 values for each sequencing run directory and store it as a json file. 
Also this script will mark the new sequencing run as ready for processing (seeded) for a specific pipeline.

## Usage
<pre><code>
find_new_seqrun_and_prepare_md5.py [-h] -p SEQRUN_PATH 
                                        -m MD5_PATH    
                                        -d DBCONFIG_PATH 
                                        -s SLACK_CONFIG 
                                        -a ASANA_CONFIG 
                                        -i ASANA_PROJECT_ID 
                                        -n PIPELINE_NAME

optional arguments:
  -h, --help              show this help message and exit
  -p /--seqrun_path       Seqrun directory path
  -m /--md5_path          Seqrun md5 output dir
  -d /--dbconfig_path     Database configuration json file
  -s /--slack_config      Slack configuration json file
  -a /--asana_config      Asana configuration json file
  -i /--asana_project_id  Asana project id
  -n /--pipeline_name     IGF pipeline name
</code></pre>

## Input file format

### DB config
<pre><code>
{"dbname":"DBNAME","driver":"sqlite"}
</code></pre>

### Asana config
<pre><code>
{ "asana_personal_token" : "XYZ" },
</code></pre>

### Slack config
<pre><code>
{"slack_token" : "ZYX", 
 "slack_channel" : "C001", 
 "slack_bot_id" : "U007"}
</code></pre>
