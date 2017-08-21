import argparse
from igf_data.process.seqrun_processing.find_and_process_new_seqrun import find_new_seqrun_dir, calculate_file_md5, load_seqrun_files_to_db

parser=argparse.ArgumentParser()
parser.add_argument('-p','--seqrun_path', required=True, help='Seqrun directory path')
parser.add_argument('-m','--md5_path', required=True, help='Seqrun md5 output dir')
parser.add_argument('-d','--dbconfig_path', required=True, help='Database configuration json file')
args=parser.parse_args()

seqrun_path=args.seqrun_path
md5_path=args.md5_path
dbconfig_path=args.dbconfig_path=

new_seqruns=find_new_seqrun_dir(seqrun_path, dbconfig_path)
new_seqrun_files_and_md5=calculate_file_md5(seqrun_info=new_seqruns, md5_out=md5_path)
load_seqrun_files_to_db(seqrun_info=new_seqruns, seqrun_md5_info=new_seqrun_files_and_md5, dbconfig=dbconfig_path)
