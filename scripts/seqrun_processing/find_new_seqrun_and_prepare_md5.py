import os, hashlib
from igf_data.igfdb.baseadaptor import BaseAdaptor
from igf_data.igfdb.seqrunadaptor import SeqrunAdaptor

def find_new_seqrun_dir(
def calculate_file_checksum(filepath, hasher='md5'):
  '''
  A method for file checksum calculation
  required param:
  filepath: a file path
  hasher: default is md5, allowed: md5 or sha256
  '''
  try:
    with open(filepath, 'rb') as infile:
      if hasher=='md5':
        file_checksum=hashlib.md5(infile.read()).hexdigest()
        return file_checksum
      elif hasher=='sha256':
        file_checksum=hashlib.sha256(infile.read()).hexdigest()
        return file_checksum
      else:
        raise('hasher {0} is not supported'.format(hasher))
  except:
    raise

if __name__=='__main__':
  checksum=calculate_file_checksum(filepath='/home/vmuser/git_code/data-management-python/doc/data/Illumina/RunInfo.xml')
  print(checksum)
