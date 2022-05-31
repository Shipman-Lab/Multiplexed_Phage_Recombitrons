from Bio import SeqIO
from determine_sequence_edits import determine_sequence_edits
import pandas as pd
import numpy as np

# manual input fastq names and search variables 

# set up fastq_dict = "fastq_name": [L_fuzzy_search, R_fuzzy_search}
#set up edits 

# are we looking for all 10 edits in every fastq? this really affects the below code
fastq_dict = {"msSBK-32-20_S20_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              "msSBK-32-43_S43_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              "msSBK-32-21_S21_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              "msSBK-32-44_S44_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              "msSBK-32-22_S22_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              "msSBK-32-45_S45_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              "msSBK-32-23_S23_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              "msSBK-32-46_S46_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              "msSBK-32-24_S24_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              "msSBK-32-52_S52_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              "msSBK-32-31_S31_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              "msSBK-32-53_S53_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              "msSBK-32-32_S32_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              "msSBK-32-55_S55_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              "msSBK-32-33_S33_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              "msSBK-32-60_S60_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              "msSBK-32-34_S34_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              "msSBK-32-61_S61_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              "msSBK-32-35_S35_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              "msSBK-32-62_S62_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              "msSBK-32-42_S42_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"]}
              #later on, can get ride of this dictionary - just define L and R fuzzy out as variables 

edits_dict = {"C4T": {"L_fuzzy_inside": "CGTCCGATAT", "R_fuzzy_inside": "ACGAANGATA"},
             "G10A": {"L_fuzzy_inside": "ATATNACGAA", "R_fuzzy_inside": "GATAAATGNA"},
             "C19T": {"L_fuzzy_inside": "ANGATAAATG", "R_fuzzy_inside": "AGCAAATGNC"},
             "C28T": {"L_fuzzy_inside": "GNAGCAAATG", "R_fuzzy_inside": "CTGAGCGGNT"},
             "T37C": {"L_fuzzy_inside": "GNCTGAGCGG", "R_fuzzy_inside": "TGTAAGTTNC"},
             "C46T": {"L_fuzzy_inside": "GNTGTAAGTT", "R_fuzzy_inside": "CGCAATAANG"},
             "C55T": {"L_fuzzy_inside": "TNCGCAATAA", "R_fuzzy_inside": "GTCGGCAANT"},
             "C64T": {"L_fuzzy_inside": "ANGTCGGCAA", "R_fuzzy_inside": "TTTGGCGGNT"},
             "C73T": {"L_fuzzy_inside": "ANTTTGGCGG", "R_fuzzy_inside": "TTCCTTTCNA"},
             "C82T": {"L_fuzzy_inside": "GNTTCCTTTC", "R_fuzzy_inside": "ATTAANAAAC"},
             "C88T": {"L_fuzzy_inside": "TTTCNATTAA", "R_fuzzy_inside": "AAACTTTCGC"}}


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
        df = pd.DataFrame(index = SeqIO.index(data_loc + fastq_name, 'fastq'),
                          columns = ["does_it_map", "correct_length"] + nt_at_edit_pos_name_array)
        L_fuzzy_search_out = fastq_dict[fastq_name][0]
        R_fuzzy_search_out = fastq_dict[fastq_name][1]
        L_fuzzy_search_in = edits_dict[edit]["L_fuzzy_inside"]
        R_fuzzy_search_in = edits_dict[edit]["R_fuzzy_inside"]
        edit_position_array = fastq_dict[fastq_name][2]
        for read in SeqIO.parse(handle, "fastq"):
            output = determine_sequence_edits(str(read.seq), L_fuzzy_search_out, R_fuzzy_search_out, L_fuzzy_search_in, R_fuzzy_search_in, edit_position_array)
            df.loc[read.id, "does_it_map"] = output['does_it_map']
            df.loc[read.id, "correct_length"] = output['correct_length']
            if output["correct_length"] == True:
                for index, col in enumerate(nt_at_edit_pos_name_array):
                    df.loc[read.id, col] = output["nt_at_edit_pos"][index]
        df.to_csv(data_loc + fastq_name[:-6] + "_processed_data.csv")


