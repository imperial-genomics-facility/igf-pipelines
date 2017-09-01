# Load platform data to db

A script for loading platform data to database

## Usage
<pre><code>load_platform_data.py [-h] -p PLATFORM_DATA [-u] -d DBCONFIG_PATH -s SLACK_CONFIG</pre></code>

## Parameters
<pre><code>  -h /--help     show this help message and exit
  -p /--platform_data    Platform data json file
  -u /--update           Update existing platform data, default: False
  -d /--dbconfig_path    Database configuration json file
  -s /--slack_config     Slack configuration json file
</code></pre>

## Platform data json file
<pre><code>{ "platform_igf_id":"ILM4K_001", 
  "model_name":"HISEQ4000",
  "vendor_name":"ILLUMINA",
  "software_name":"RTA",
  "software_version":"RTA2"}</code></pre>
