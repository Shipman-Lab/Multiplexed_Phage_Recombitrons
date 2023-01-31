"""
BEFORE BEGINNING: mount the Shipman-Lab hive drive on your Mac (what is currently supported by this script)

TO RUN:
    change everything between the two hashes to match your current run
        (note: you must include your run information in the file_key.xlsx file tracked by GitHub)
    commit this file after saving & changes - the code will not run with uncommitted changes
    navigate to this folder in terminal
    python3 -m edit_site_analysis_script
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
import os
import numpy as np
from edit_site_analysis_functions import extract_and_match

def run_single_nt_edit_analysis(run_path, run_name, file_key_path, fuzziness):
    if not os.path.isdir(run_path):
        raise ValueError("Input run path is not a directory!")

    # check git hash & that there are no uncommitted changes
    status = subprocess.check_output(["git", "status"])
    if "Changes not staged for commit" in str(status, 'utf-8').strip():
        raise ValueError("Uncommitted changes - please commit before running")
    git_short_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])
    git_short_hash = str(git_short_hash, "utf-8").strip()

    # load in file key
    # if you don't care about some of these, just leave them blank in the file key
    file_key = pd.read_excel(file_key_path, converters={'barcode_num': str})
    file_key = file_key.set_index("barcode_num")

    outcome_df = file_key[["construct", "rep", "type", "gene",
                           "plasmid", "direction", "edit_name",
                           "genome_position"]].copy()
    outcome_df[["total_num_reads", "wt", "edited", "unmatched_region", "unmatched_edit_nt"]] = np.NaN

    for root, dirs, files in os.walk(run_path):
        # find folder with run_id
        for directory in dirs:
            for barcode in file_key.index:
                construct = file_key.loc[barcode, "construct"]
                rep = file_key.loc[barcode, "rep"]
                if barcode in directory:
                    print("----------------------------")
                    print("working on construct %s, rep %s"%(construct, rep))
                    start_time = time.time()
                    barcode_dir = os.path.join(root, directory)
                    fastqs = os.listdir(barcode_dir)

                    L_inside = file_key.loc[barcode, "L_inside"]
                    R_inside = file_key.loc[barcode, "R_inside"]
                    L_outside = file_key.loc[barcode, "L_outside"]
                    R_outside = file_key.loc[barcode, "R_outside"]
                    wt_nt = file_key.loc[barcode, "wt_nt"]
                    edited_nt = file_key.loc[barcode, "edited_nt"]

                    outcomes_dict = {'wt':0, 'edited':0, 'unmatched_region':0, 'unmatched_edit_nt':0}
                    for fastq in fastqs:
                        if ("fastq.gz" not in fastq):
                            continue
                        full_path = os.path.join(barcode_dir, fastq)
                        with gzip.open(full_path, 'rb') as f_in:
                            with open(full_path[:-3], 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                        records = list(SeqIO.parse(full_path[:-3], "fastq"))

                        for read in records:
                            outcomes_dict[extract_and_match(read.seq, L_outside, R_outside, L_inside,
                                                            R_inside, wt_nt, edited_nt, fuzziness=fuzziness)] += 1
                    # put into output df
                    outcome_df.loc[barcode, "total_num_reads"] = len(records)
                    outcome_df.loc[barcode, ["wt", "edited", "unmatched_region", "unmatched_edit_nt"]] = outcomes_dict
                    outcome_df.to_excel("%s_summary_df_vers"%(run_name) + str(git_short_hash) + ".xlsx")
                    print("---  processing took %s seconds ---" % (time.time() - start_time))
