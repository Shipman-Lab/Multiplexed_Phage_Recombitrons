"""
BEFORE BEGINNING: mount the Shipman-Lab hive drive on your Mac (what is currently supported by this script)
"""

##Import
from Bio import SeqIO
import pandas as pd
from collections import Counter
import time
import subprocess
import gzip
import shutil
import time

# check git hash & that there are no uncommitted changes
status = subprocess.check_output(["git", "status"])
if "Changes not staged for commit" in str(status, 'utf-8').strip():
    raise ValueError("Uncommitted changes - please commit before running")
git_short_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])
git_short_hash = str(git_short_hash, "utf-8").strip()

if not os.path.isdir("/Volumes/Shipman-Lab/BaseSpace"):
    raise ValueError("Make sure to mount the Shipman-Lab hive drive")

# load in file key
# if you don't care about some of these, just leave them blank in the file key
file_key = pd.read_excel("file_key.xlsx")\
           [["phage", "gene", "plasmid", "direction",
           "edit_name", "genome_position", "wt_nt", 
           "edited_nt", "L_inside", "R_inside", 
           "L_outside", "R_outside",
           "rep_1", "rep_2", "rep_3",
           "rep_4", "rep_5"]]
run_path = "/Volumes/Shipman-Lab/BaseSpace/msKDC_01-353873095/FASTQ_Generation_2022-06-04_18_18_55Z-570470901"
outcome_df = file_key.melt(id_vars=["phage", "gene", "plasmid", "direction",
                                    "edit_name", "genome_position", "wt_nt",
                                    "edited_nt", "L_inside", "R_inside", "L_outside", "R_outside"],
                           value_vars=["rep_1", "rep_2", "rep_3", "rep_4", "rep_5"])
# assumes you have 9 or fewer replicates per type!!!
outcome_df["rep"] = outcome_df["variable"].str[-1]
outcome_df = outcome_df.rename(columns={"value": "run_id"})
outcome_df[["wt", "edited", "unmatched_region", "unmatched_edit_nt"]] = np.NaN

for root, dirs, files in os.walk(run_path):
    for directory in dirs:
        files = os.listdir(os.path.join(root, directory))
        for name in files:
            if "fastq.gz" not in name:
                continue
            # match to file key
            for fastq_name in outcome_df["run_id"]:
                if fastq_name in name:
                    print("working on %s" %fastq_name)
                    start_time = time.time()
                    outcomes_dict = {'wt':0,'edited':0,'unmatched_region':0,'unmatched_edit_nt':0}
                    all_reads_str = []
                    read_counter = []
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
                                ### NEED TO FIGURE OUT BELOW LINE!
                                import pdb
                                pdb.set_trace()
                                outcomes_dict[extract_and_match(read, i, rep)] += read_counter[read]
                        # put into output df
                print("---  processing took %s seconds ---" % (time.time() - start_time))
                outcome_df.to_excel("msKDC001_summary_df_vers" + str(git_short_hash) + ".xlsx")


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