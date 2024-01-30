import pandas as pd
from Bio import SeqIO
import gzip
import os, sys

def demultiplex_to_list(run,file,rep,rep_bc):
    file_folder, file_path = get_file_path(run,file)
    return parse_fastq_gz(file_path,rep_bc)

def demultiplex_to_file(run,file,rep,rep_bc):
    file_folder, file_path = get_file_path(run,file)
    seqs = parse_fastq_gz(file_path,rep_bc)
    SeqIO.write(seqs,'%s\%s_%s.fastq' % (file_folder,file,rep), "fastq")

def parse_fastq_gz(filename,bc):
    seq_list = []
    with gzip.open(filename, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if bc in str(record.seq[:43]):
                seq_list.append(record)
    return seq_list
            
def get_file_path(run,file):
    miseq_folder_names = os.listdir('Y:\BaseSpace')
    miseq_folder_dict = {}
    for folder in miseq_folder_names:
        miseq_folder_dict[folder.split('-')[0]] = folder
    subfolder_names = os.listdir('Y:\BaseSpace\%s' % miseq_folder_dict[run])
    sub_folder_dict = {}
    for folder in subfolder_names:
        try: sub_folder_dict[folder.split('_')[2]] = folder
        except IndexError: continue
    file_folder = 'Y:\BaseSpace\%s\%s' % (miseq_folder_dict[run],sub_folder_dict[file.split('_')[2]])
    file_names = os.listdir(file_folder)
    name_dict = {}
    for name in file_names:
        try: name_dict[name.split('_')[5]] = name  #just dealing with forward reads
        except IndexError: continue
    file_path = 'Y:\BaseSpace\%s\%s\%s' % (miseq_folder_dict[run],sub_folder_dict[file.split('_')[2]],name_dict['R1'])
    return file_folder, file_path