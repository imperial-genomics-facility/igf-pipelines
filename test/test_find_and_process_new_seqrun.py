import unittest, json
from sqlalchemy import create_engine
from igf_data.igfdb.baseadaptor import BaseAdaptor
from igf_data.igfdb.platformadaptor import PlatformAdaptor
from igf_data.igfdb.seqrunadaptor import SeqrunAdaptor
from igflib.process.seqrun_processing.find_and_process_new_seqrun import find_new_seqrun_dir,check_finished_seqrun_dir,check_seqrun_dir_in_db,calculate_file_md5,calculate_file_checksum,prepare_seqrun_for_db,load_seqrun_files_to_db

class Find_seqrun_test1(unittest.TestCase):
  def setup(self):
    self.path='../data/seqrun_dir'
    self.dbconfig='../data/dbconfig.json'
    self.md5_out_path='../data'
    seqrun_json='../data/seqrun_db_data.json'
    platform_json='../data/platform_db_data.json'

    dbparam=None
    with open(dbconfig, 'r') as json_data:
      dbparam=json.load(json_data)
    base=BaseAdaptor(**dbparam)
    self.engine=base.engine
    Base.metadata.create_all(self.engine)
    base.start_session()

    with open(platform_json, 'r') as json_data:
      platform_data=json.load(json_data)
      pl=PlatformAdaptor(**{'session':base.session})
      pl.store_platform_data(data=platform_data)
        
   with open(seqrun_json, 'r') as json_data:
     seqrun_data=json.load(json_data)
     sra=SeqrunAdaptor(**{'session':base.session})
     sra.store_seqrun_and_attribute_data(data=seqrun_data)
    base.stop_session()

  def tearDown(self):
     Base.metadata.create_all(self.engine)

   
