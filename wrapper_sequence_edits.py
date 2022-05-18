from Bio import SeqIO
from determine_sequence_edits import determine_sequence_edits
import pandas as pd
import numpy as np

# manual input fastq names and search variables 

# set up dictionary = "fastq_name": [L_fuzzy_search, R_fuzzy_search, L_to_R_length, edit_position_array]}
# edit_position_array = [base_position1, base_position2, etc.] or [4, 15, 27] which are all referred to as base_position in determine_sequencing_edits

fastq_dict = {"msSBK-29-37_S37_L001_R1_001.fastq": ["TGTAAGTTCCGCAATAACGT", "CGATGTGCGCCAGCGGAGTC", 90, [45, 62]]}
             

# direct it to where the reads are
data_loc = "./data/"  

# handle is the opened fastq file
# read is each fastq in the iterator 
# iterator (big list that you can loop through) is made by SeqIO.parse 

# open fastq of interest
# define parameters of dataframe to create 
# send each read through determine_sequence_edits.py 
# put into output csv from the dataframe and loop back to next read 

for fastq_name in fastq_dict.keys():
    print(fastq_name)
    with open(data_loc+fastq_name) as handle:
        nt_at_edit_pos_name_array = ["nt_at_edit_pos_%i"%(x) for x in fastq_dict[fastq_name][3]]
        df = pd.DataFrame(index = SeqIO.index(data_loc + fastq_name, 'fastq'),
                          columns = ["does_it_map", "correct_length"] + nt_at_edit_pos_name_array)
        L_fuzzy_search = fastq_dict[fastq_name][0]
        R_fuzzy_search = fastq_dict[fastq_name][1]
        L_to_R_length = fastq_dict[fastq_name][2]
        edit_position_array = fastq_dict[fastq_name][3]
        for read in SeqIO.parse(handle, "fastq"):
            output = determine_sequence_edits(str(read.seq), L_fuzzy_search, R_fuzzy_search, L_to_R_length, edit_position_array)
            df.loc[read.id, "does_it_map"] = output['does_it_map']
            df.loc[read.id, "correct_length"] = output['correct_length']
            if output["correct_length"] == True:
                for index, col in enumerate(nt_at_edit_pos_name_array):
                    df.loc[read.id, col] = output["nt_at_edit_pos"][index]
        df.to_csv(data_loc + fastq_name[:-6] + "_processed_data.csv")


