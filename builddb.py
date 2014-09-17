#!/usr/bin/env python
import subprocess
from os import listdir 
from os.path import join
import re

prefix = './'
file_regex = re.compile('.fa$')
file_list = [f for f in listdir(prefix) if (file_regex.search(f) != None)]
for fasta_file in file_list:
  if prefix == './':
    fasta_path = fasta_file
  else:
    fasta_path = join(prefix, fasta_file);
  species = fasta_file[:fasta_file.find('.')]
  print 'Building db for ', species
  subprocess.call(['formatdb', '-i', fasta_path, '-o', 'T', '-n', species])

