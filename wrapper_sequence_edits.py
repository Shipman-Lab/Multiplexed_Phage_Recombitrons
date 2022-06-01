from Bio import SeqIO
from determine_sequence_edits import determine_sequence_edits
import pandas as pd
import numpy as np 

# import pdb
# pdb.set_trace()

# manual input fastq names and which edits to look for 

#set up edits 
#this is for edits across the internal L-tail editor on Lambda (L5, L4, L3, L2, L1, center, R1, R2, R3, R4, R5)
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
             #since some fuzzies have N at the adjacent edit site, a fuzziness of 1 is necessary to enable a match and 2 would be more lenient 

#recode parameters 
central_edits_dict = {key: edits_dict[key] for key in ("C46T",)}
L1_R1_edits_dict = {key: edits_dict[key] for key in ("T37C", "C46T", "C55T")}
L2_R2_edits_dict = {key: edits_dict[key] for key in ("C28T", "T37C", "C46T", "C55T", "C64T")}
L3_R3_edits_dict = {key: edits_dict[key] for key in ("C19T", "C28T", "T37C", "C46T", "C55T", "C64T", "C73T")}
L4_R4_edits_dict = {key: edits_dict[key] for key in ("G10A", "C19T", "C28T", "T37C", "C46T", "C55T", "C64T", "C73T", "C82T")}
L5_R5_edits_dict = {key: edits_dict[key] for key in ("C4T", "G10A", "C19T", "C28T", "T37C", "C46T", "C55T", "C64T", "C73T", "C82T", "C88T")}

#for homology, MOI, and temperature parameters use central_edits_dict (pSBK.164)

#for edit placement, add keys individually 

#outside fuzzy search for internal L-tail editor on lambda 
L_fuzzy_search_out = "TATGACCAGCCAACGTCCGA"
R_fuzzy_search_out = "ACTTTCGCAGTAAATCCCAT"


#set up fastq_dict = {"fastq_name": [L_fuzzy_search_out, R_fuzzy_search_out]}
fastq_dict = {"msSBK-32-20_S20_L001_R1_001_fortest.fastq": [L_fuzzy_search_out, R_fuzzy_search_out, central_edits_dict]}
              # "msSBK-32-43_S43_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              # "msSBK-32-21_S21_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              # "msSBK-32-44_S44_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              # "msSBK-32-22_S22_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              # "msSBK-32-45_S45_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              # "msSBK-32-23_S23_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              # "msSBK-32-46_S46_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              # "msSBK-32-24_S24_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              # "msSBK-32-52_S52_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              # "msSBK-32-31_S31_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              # "msSBK-32-53_S53_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              # "msSBK-32-32_S32_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              # "msSBK-32-55_S55_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              # "msSBK-32-33_S33_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              # "msSBK-32-60_S60_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              # "msSBK-32-34_S34_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              # "msSBK-32-61_S61_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              # "msSBK-32-35_S35_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],   
              # "msSBK-32-62_S62_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"],
              # "msSBK-32-42_S42_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT"]}

# #IF ANALYZING CI MUT EDITING 
# #outside fuzzy search for cI mut editor on lambda 
# L_cI_out = "GGCAGTCAGGCGTTGGTGCT"
# R_cI_out = "ATCGCCAGAGAAATCTACGA"

# edits_dict_cI = {"T23A": {"L_fuzzy_inside": "ATCAATGCAT", "R_fuzzy_inside": "AAATGCTTAT"},
#                  "A52T": {"L_fuzzy_inside": "ATTGCTTGCA", "R_fuzzy_inside": "AAATTCTCAA"}}


# direct it to where the reads are
data_loc = "./data/"  

# handle is the opened fastq file
# index is all of the reads that are in the fastq 
# read is each fastq in the iterator 
# iterator (big list that you can loop through) is made by SeqIO.parse 

# open fastq of interest
# define parameters of dataframe to create 
# send each read through determine_sequence_edits.py 
# put into output csv from the dataframe and loop back to next read 


for fastq_name in fastq_dict.keys():
    print(fastq_name)
    with open(data_loc+fastq_name) as handle:
        fastq_edits_dict = fastq_dict[fastq_name][2]
        nt_col_names = ["%s_nt"%edit_name for edit_name in fastq_edits_dict.keys()]
        map_col_names = ["%s_map"%edit_name for edit_name in fastq_edits_dict.keys()]
        df = pd.DataFrame(index = SeqIO.index(data_loc + fastq_name, 'fastq'),
                          columns = ["does_it_map"] + nt_col_names + map_col_names)
        L_fuzzy_search_out = fastq_dict[fastq_name][0]
        R_fuzzy_search_out = fastq_dict[fastq_name][1]
        for read in SeqIO.parse(handle, "fastq"):
            output = determine_sequence_edits(str(read.seq), L_fuzzy_search_out, R_fuzzy_search_out, fastq_edits_dict)
            df.loc[read.id, "does_it_map"] = output['does_it_map']
            does_it_map_around_edit_dict = output['does_it_map_around_edit_dict']
            nt_at_edit_pos_dict = output['nt_at_edit_pos_dict']
            for edit_name in does_it_map_around_edit_dict.keys():
                df.loc[read.id, edit_name + "_nt"] = nt_at_edit_pos_dict[edit_name]
                df.loc[read.id, edit_name + "_map"] = does_it_map_around_edit_dict[edit_name]
        df.to_csv(data_loc + fastq_name[:-6] + "_processed_data.csv")


