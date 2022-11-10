##Import
from Bio import SeqIO
import pandas as pd
from collections import Counter
import time
import subprocess
import gzip
import shutil
import time
import os
import argparse
import numpy as np
from edit_site_analysis_functions import extract_and_match

# check git hash & that there are no uncommitted changes
status = subprocess.check_output(["git", "status"])
if "Changes not staged for commit" in str(status, 'utf-8').strip():
    raise ValueError("Uncommitted changes - please commit before running")
git_short_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])
git_short_hash = str(git_short_hash, "utf-8").strip()

parser = argparse.ArgumentParser()
parser.add_argument("input", help = "Path to input excel file")
parser.add_argument("fastq", help = "Path to directory containing sequencing reads")
parser.add_argument("output", help = "Path to directory containing output csv")
args = parser.parse_args()

# load in edits key
# if you don't care about some of these, just leave them blank in the edits key
edits_key = pd.read_excel("edits_key.xlsx")\
           [["edit_name", "strain", "gene", "plasmid",
           "direction", "genome_position", "wt_nt", 
           "edited_nt", "L_inside", "R_inside", 
           "L_outside", "R_outside"]]
# load in file key
# multiple edit names are grouped in the same cell, delimited by commas
file_key = pd.read_excel(args.input)[["run_id", "edit_name"]]
# create output list (each row of the list will be results for a unique run_id + edit_name pair)
output_list = []

for root, dirs, files in os.walk(args.fastq):
    for directory in dirs:
        files = os.listdir(os.path.join(root, directory))
        for name in files:
            if "fastq.gz" not in name:
                continue
            # match to file key
            for row in file_key.iterrows():
                run_id = row[1]['run_id']
                if run_id in name or '-'.join(run_id.split('_')) in name:
                    print("working on %s" % run_id)
                    start_time = time.time()
                    for edit_name in [x.strip() for x in row[1]['edit_name'].split(',')]: 
                        edit_ix = edits_key[edits_key['edit_name'] == edit_name].index[0]
                        outcomes_dict = {'wt':0, 'edited':0, 'unmatched_region':0, 'unmatched_edit_nt':0}
                        all_reads_str = []
                        read_counter = []
                        L_outside = edits_key.iloc[edit_ix]['L_outside']
                        R_outside = edits_key.iloc[edit_ix]['R_outside']
                        L_inside = edits_key.iloc[edit_ix]['L_inside']
                        R_inside = edits_key.iloc[edit_ix]['R_inside']
                        wt_nt = edits_key.iloc[edit_ix]['wt_nt']
                        edited_nt = edits_key.iloc[edit_ix]['edited_nt']
                        # need to go into the folder & unzip the file
                        full_path = os.path.join(root, directory, name)
                        with gzip.open(full_path, 'rb') as f_in:
                            with open(full_path[:-3], 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                            with open(full_path[:-3]) as handle:
                                for record in SeqIO.parse(handle, "fastq"):
                                    all_reads_str.append(str(record.seq))
                                read_counter = Counter(all_reads_str)
                                for read in read_counter:
                                    outcomes_dict[extract_and_match(read, L_outside, R_outside, L_inside,
                                                                    R_inside, wt_nt, edited_nt)] += read_counter[read]
                        # put into output list
                        output_list.append([run_id, edit_name, edits_key.iloc[edit_ix]['strain'], 
                                            edits_key.iloc[edit_ix]['gene'], edits_key.iloc[edit_ix]['plasmid'],
                                            edits_key.iloc[edit_ix]['direction'], edits_key.iloc[edit_ix]['genome_position'], 
                                            wt_nt, edited_nt, L_inside, R_inside, L_outside, R_outside, len(all_reads_str),
                                            outcomes_dict['wt'], outcomes_dict['edited'], outcomes_dict['unmatched_region'],
                                            outcomes_dict['unmatched_edit_nt'], 100.*outcomes_dict['edited']/len(all_reads_str)])
                    print("---  processing took %s seconds ---" % (time.time() - start_time))

# create output df
output_df = pd.DataFrame(output_list, columns=["run_id", "edit_name", "strain", "gene", "plasmid",
                                    "direction", "genome_position", "wt_nt",
                                    "edited_nt", "L_inside", "R_inside", "L_outside", "R_outside", "total_reads",
                                    "wt", "edited", "unmatched_region", "unmatched_edit_nt", "perc_edited"])
#output_df.to_csv(os.path.join(args.output, "test_output.csv"))
outcome_df.to_csv(os.path.join(args.output, "summary_df_vers" + str(git_short_hash) + ".csv"))


# for i in samples.index:
#     sample_i = int(i)
#     outcomes_dict[i] = {}
#     for rep in ["rep_1", "rep_2", "rep_3", "rep_4", "rep_5"]:
#         outcomes_dict[i][rep]= {'wt':0,'edited':0,'unmatched_region':0,'unmatched_edit_nt':0}
#         all_reads_str = []
#         read_counter = []
#         fastq_reads = "./%s_trimmed.fastq" % samples[rep][i]
#         try:
#             for seq_record in SeqIO.parse(fastq_reads, "fastq"):
#                 all_reads_str.append(str(seq_record.seq))
#             read_counter = Counter(all_reads_str)
#             for read in read_counter:
#                 outcomes_dict[sample_i][rep][extract_and_match(read,i,rep)] += read_counter[read]
#             print(samples[rep][i])
#         except IOError: #this happens when a file is missing
#             print("%s missing" % samples[rep][i])