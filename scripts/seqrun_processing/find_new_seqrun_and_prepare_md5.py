import argparse
from igf_data.task_tracking.igf_slack import IGF_slack
from igf_data.task_tracking.igf_asana import IGF_asana
from igf_data.process.seqrun_processing.find_and_process_new_seqrun import find_new_seqrun_dir, calculate_file_md5, load_seqrun_files_to_db, seed_pipeline_table_for_new_seqrun

parser=argparse.ArgumentParser()
parser.add_argument('-p','--seqrun_path', required=True, help='Seqrun directory path')
parser.add_argument('-m','--md5_path', required=True, help='Seqrun md5 output dir')
parser.add_argument('-d','--dbconfig_path', required=True, help='Database configuration json file')
parser.add_argument('-s','--slack_config', required=True, help='Slack configuration json file')
parser.add_argument('-a','--asana_config', required=True, help='Asana configuration json file')
parser.add_argument('-i','--asana_project_id', required=True, help='Asana project id')
parser.add_argument('-n','--pipeline_name', required=True, help='IGF pipeline name')
args=parser.parse_args()

seqrun_path=args.seqrun_path
md5_path=args.md5_path
dbconfig_path=args.dbconfig_path
slack_config=args.slack_config
asana_config=args.asana_config
asana_project_id=args.asana_project_id
pipeline_name=args.pipeline_name

slack_obj=IGF_slack(slack_config=slack_config)
asana_obj=IGF_asana(asana_config=asana_config, asana_project_id=asana_project_id)

try:
  new_seqruns=find_new_seqrun_dir(seqrun_path, dbconfig_path)
  if len(new_seqruns.keys()) > 0:
    new_seqrun_files_and_md5=calculate_file_md5(seqrun_info=new_seqruns, md5_out=md5_path,seqrun_path=seqrun_path)
    load_seqrun_files_to_db(seqrun_info=new_seqruns, seqrun_md5_info=new_seqrun_files_and_md5, dbconfig=dbconfig_path)
    seed_pipeline_table_for_new_seqrun(seqrun_info=new_seqruns, pipeline_name=pipeline_name, dbconfig=dbconfig_path)

    for seqrun_name in new_seqruns.keys():
      message='found new sequencing run {0}'.format(seqrun_name)
      res=asana_obj.comment_asana_task(task_name=seqrun_name, comment=message)
      slack_obj.post_message_to_channel(message,reaction='pass')
      message='New asana task created for seqrun {0}, url: https://app.asana.com/0/{1}/{2}'.format(seqrun_name, asana_project_id, res['target']['id'])
      slack_obj.post_message_to_channel(message,reaction='pass')
  else:
    slack_obj.post_message_to_channel(message='No new sequencing run found',reaction='sleep')
except Exception as e:
  message='Failed to load new seqruns, received following error: {0}'.format(e)
  slack_obj.post_message_to_channel(message,reaction='fail')
  raise
