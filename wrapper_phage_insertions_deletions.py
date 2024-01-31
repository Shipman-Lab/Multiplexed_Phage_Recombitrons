# this is a script to process all fastqs in a subfolder named "data"
# it will save two outputs per fastq
#
# one output is fastq_name_map_dict.p
#   this pickle contains a dictionary with four numbers:
#       total number of reads processed
#       total number of reads mapped (based on the fuzzy search parameters)
#       total number of reads containing a barcode
#       total reads with code errors (currently only containing 2 barcodes)
#
# the other output is a csv of a dataframe
#   the index is the barcode
#   the value is the count of that barcode identified in the fastq

import pandas as pd
import pickle as p
import numpy as np
import os
import gzip
import shutil
import time
import subprocess
from Bio import SeqIO
from phage_insertions_deletions import amplicon_insertion_deletion, shotgun_insertion_deletion

# check git hash & that there are no uncommitted changes
status = subprocess.check_output(["git", "status"])
if "Changes not staged for commit" in str(status, 'utf-8').strip():
    raise ValueError("Uncommitted changes - please commit before running")
git_short_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])
git_short_hash = str(git_short_hash, "utf-8").strip()

run_path = os.path.expanduser("~/Desktop/msAGK_08-392858067/FASTQ_Generation_2023-07-11_10_23_24Z-681323663/")
if not os.path.isdir(run_path):
    raise ValueError("Input run path is not a real directory!")

file_key = pd.read_excel("editing_file_key.xlsx", index_col="file_id")

# left & right flank should be out of any deletion regions
right_flank = "CATACATCTCGTAGATTTCT"

summary_df = pd.DataFrame(index=file_key.index, columns=["type", "indel_len", "num_reads",
                                                         "fraction_mapped",
                                                         "fraction_w_middle_search_seq",
                                                         "fraction_middle_search_error"])

for root, dirs, files in os.walk(run_path):
    for directory in dirs:
        files = os.listdir(os.path.join(root, directory))
        for name in files:
            if "fastq.gz" not in name:
                continue
            # match to file key
            for fastq_name in file_key.index:
                if fastq_name in name:
                    print("working on %s" %fastq_name)
                    summary_df.loc[fastq_name, "type"] = file_key.loc[fastq_name, "type"]
                    summary_df.loc[fastq_name, "indel_len"] = len(file_key.loc[fastq_name,"middle_search_seq"])
                    start_time = time.time()
                    # need to go into the folder & unzip the file
                    full_path = os.path.join(root, directory, name)
                    with gzip.open(full_path, 'rb') as f_in:
                        with open(full_path[:-3], 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                        output_df = pd.DataFrame(0, index= SeqIO.index(full_path[:-3], "fastq"),
                                                 columns=["does_it_map", "contain_middle_search"])
                        with open(full_path[:-3]) as handle:
                            for record in SeqIO.parse(handle, "fastq"):
                                seq = record.seq 
                                read_output = amplicon_insertion_deletion(seq, file_key.loc[fastq_name, "left_flank"], right_flank,
                                                                          file_key.loc[fastq_name, "middle_search_seq"],
                                                                          file_key.loc[fastq_name, "fuzziness"])
                                output_df.loc[record.id, "does_it_map"] = read_output["does_it_map"]
                                output_df.loc[record.id, "contain_middle_search"] = read_output["contain_middle_search"]
                    summary_df.loc[fastq_name, "num_reads"] = len(output_df["does_it_map"])
                    summary_df.loc[fastq_name, "fraction_mapped"] = sum(output_df["does_it_map"])/len(output_df["does_it_map"])
                    summary_df.loc[fastq_name, "fraction_middle_search_error"] = sum(output_df["contain_middle_search"] == "Error")/len(output_df["contain_middle_search"])
                    summary_df.loc[fastq_name, "fraction_w_middle_search_seq"] = sum(output_df["contain_middle_search"] == True)/len(output_df["contain_middle_search"])
                    output_df.to_excel("./processed/" + fastq_name + "_read_counts_vers" + str(git_short_hash) + ".xlsx")
                    print("---  processing took %s seconds ---" % (time.time() - start_time))
                    summary_df.to_excel("msAGK08_summary_df.xlsx")

# nextera_dict = {"msSBK-32-205_S205_L001_R1_001.fastq": [delete_4bp, 0],
#                  "msSBK-32-206_S206_L001_R1_001.fastq": [delete_4bp, 0],
#                  "msSBK-32-207_S207_L001_R1_001.fastq": [delete_4bp, 0],
#                  "msSBK-32-208_S208_L001_R1_001.fastq": [delete_32bp, 5],
#                  "msSBK-32-209_S209_L001_R1_001.fastq": [delete_32bp, 5],
#                  "msSBK-32-210_S210_L001_R1_001.fastq": [delete_32bp, 5],
#                  "msSBK-32-211_S211_L001_R1_001.fastq": [insert_4bp, 0],
#                  "msSBK-32-212_S212_L001_R1_001.fastq": [insert_4bp, 0],
#                  "msSBK-32-213_S213_L001_R1_001.fastq": [insert_4bp, 0],
#                  "msSBK-32-214_S214_L001_R1_001.fastq": [insert_34bp, 5],
#                  "msSBK-32-215_S215_L001_R1_001.fastq": [insert_34bp, 5],
#                  "msSBK-32-216_S216_L001_R1_001.fastq": [insert_34bp, 5]}

# for fastq_name in nextera_dict.keys():
#     print("working on %s" %fastq_name)
#     start_time = time.time()
#     with open("./data/" + fastq_name) as handle:
#         output_df = pd.DataFrame(0, index=SeqIO.index("./data/" + fastq_name, "fastq"),
#                                  columns=["does_it_map", "contain_middle_search"])
#         for record in SeqIO.parse(handle, "fastq"):
#             seq = record.seq 
#             read_output = shotgun_insertion_deletion(seq, left_flank, right_flank,
#                                                       nextera_dict[fastq_name][0],
#                                                       nextera_dict[fastq_name][1])
#             output_df.loc[record.id, "does_it_map"] = read_output["does_it_map"]
#             output_df.loc[record.id, "contain_middle_search"] = read_output["contain_middle_search"]

#     output_df.to_excel(fastq_name + "_read_counts.xlsx")
#     print("---  processing took %s seconds ---" % (time.time() - start_time))


