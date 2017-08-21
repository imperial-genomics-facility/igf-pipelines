from igflib.process.seqrun_processing.find_and_process_new_seqrun import find_new_seqrun_dir, calculate_file_md5, load_seqrun_files_to_db

if __name__=='__main__':
  path='/home/vmuser/git_code/igf-pipelines/data/seqrun_dir'
  md5_out_path='/home/vmuser/git_code/test_dir'
  #checksum=calculate_file_checksum(filepath='/home/vmuser/git_code/data-management-python/doc/data/Illumina/RunInfo.xml')
  #print(checksum)

  dbconfig='/home/vmuser/git_code/igf-pipelines/data/dbconfig.json'
  new_seqrun=find_new_seqrun_dir(path, dbconfig)
  new_seqrun_and_md5=calculate_file_md5(seqrun_info=new_seqrun, md5_out=md5_out_path)
  load_seqrun_files_to_db(seqrun_info=new_seqrun, seqrun_md5_info=new_seqrun_and_md5, dbconfig=dbconfig)
  

