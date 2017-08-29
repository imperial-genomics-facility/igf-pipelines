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
  -h, --help            show this help message and exit
  -p /--seqrun_path      SEQRUN_PATH     :  Seqrun directory path
  -m /--md5_path         MD5_PATH        :  Seqrun md5 output dir
  -d /--dbconfig_path    DBCONFIG_PATH   :  Database configuration json file
  -s /--slack_config     SLACK_CONFIG    :  Slack configuration json file
  -a /--asana_config     ASANA_CONFIG    :  Asana configuration json file
  -i /--asana_project_id ASANA_PROJECT_ID:  Asana project id
  -n /--pipeline_name    PIPELINE_NAME   :  IGF pipeline name

</code></pre>

## Input file format

### DB config
<code>
{"dbname":"DBNAME","driver":"sqlite"}
</code>

### Asana config
<code>
{ "asana_personal_token" : "XYZ" },
</code>

### Slack config
<code>
{"slack_token" : "ZYX", 
 "slack_channel" : "C001", 
 "slack_bot_id" : "U007"}
</code>
