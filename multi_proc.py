#!/usr/bin/env python
#$ -cwd
#$ -N blast_4th
#$ -pe shm 6
import multiprocessing
from multiprocessing import Process
import subprocess
import os
from os import listdir 
from os.path import join
from decimal import *
from time import time

def blastwork(species, prefix, file_list, out_fname, num_threads=3):
    with open(out_fname, 'w') as out_file:
        for f in file_list:
            fname = join(prefix, f)
            blast_pipe = subprocess.Popen(['blastp', '-db', \
                                           species, '-query', \
                                           fname, '-outfmt', '6', \
                                           '-max_target_seqs', '1', \
                                           '-num_threads', str(num_threads)], 
                                           stdout=subprocess.PIPE)
            (stdoutdata, stdindata) = blast_pipe.communicate()
            blast_data = stdoutdata.split('\t')
            if len(blast_data) > 1:
                blast_score = blast_data[11]
            else:
                blast_score = '0\n'
            out_file.write('{}\t{}'.format(f, blast_score))

def blastdb(species, prefix, num_proc=6):
    file_list = [f for f in listdir(prefix) if f.endswith('.fa') ]
    proc_work = int(len(file_list) / 6)
    start = 0
    proc_list = []
    for i in range(num_proc):
        end = min(len(file_list), (i + 1) * proc_work)
        sub_list = file_list[start:end]
        start = end
        db_info = os.stat(species + '.phr')
#        num_thread = int(math.sqrt(db_info.st_size / 600000.0))
#        num_thread = max(num_thread, 3)
        num_thread = 5
        out_fname = species + str(i)
        proc_list.append(Process(target=blastwork, args=(species, prefix, sub_list, out_fname,num_thread,)))
        proc_list[i].start()
    for p in proc_list:
        p.join()
    out_fname = species + '.blastp'
    with open(out_fname, 'w') as out_file:
        for i in range(num_proc):
            in_fname = species + str(i)
            with open(in_fname) as in_file:
                for line in in_file:
                    out_file.write(line)
            os.remove(in_fname)

if __name__ == '__main__':
    prefix = './'
    species_list = [f for f in listdir(prefix) if f.endswith('.phr') ]
    t0 = time()
    for spc in species_list:
        species = spc[:-4]
        print 'BLAST against {}...'.format(species)
        if os.path.isfile(species + '.done'):
            print 'Already done.\nSkipping...'
            continue
        if os.path.isfile(species + '0'):
            print 'Other instance working on it.\nSkipping...'
            continue
        with open(species + '.done', 'w') as rec_file:
            rec_file.write('Working on this\n')
        
        blastdb(species, 'down/')
        with open(species + '.done', 'w') as rec_file:
            rec_file.write('{}\n'.format(species))
        t1 = time()
        print 'Time elapsed {}s'.format(t1 - t0)
        if t1 - t0 >= 40 * 3600:
            break
    print 'End process...'
