from Bio import SeqIO
from determine_sequence_edits import determine_sequence_edits
import pandas as pd
import numpy as np 

# import pdb
# pdb.set_trace()

# manual input fastq names and which edits to look for 

#set up edits 
#this is for edits across the internal L-tail editor on Lambda (L5, L4, L3, L2, L1, center, R1, R2, R3, R4, R5) + A45G is for the pSBK.160 dRT control
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
             "C88T": {"L_fuzzy_inside": "TTTCNATTAA", "R_fuzzy_inside": "AAACTTTCGC"},
             "A45G": {"L_fuzzy_inside": "CTTTCGCAGT", "R_fuzzy_inside": "AATCCCATGA"}}
             #since some fuzzies have N at the adjacent edit site, a fuzziness of 1 is necessary to enable a match and 2 would be more lenient 

#recode parameters
central_edits_dict = {key: edits_dict[key] for key in ("C46T",)}
L1_R1_edits_dict = {key: edits_dict[key] for key in ("T37C", "C46T", "C55T")}
L2_R2_edits_dict = {key: edits_dict[key] for key in ("C28T", "T37C", "C46T", "C55T", "C64T")}
L3_R3_edits_dict = {key: edits_dict[key] for key in ("C19T", "C28T", "T37C", "C46T", "C55T", "C64T", "C73T")}
L4_R4_edits_dict = {key: edits_dict[key] for key in ("G10A", "C19T", "C28T", "T37C", "C46T", "C55T", "C64T", "C73T", "C82T")}
L5_R5_edits_dict = {key: edits_dict[key] for key in ("C4T", "G10A", "C19T", "C28T", "T37C", "C46T", "C55T", "C64T", "C73T", "C82T", "C88T")}
dRT_edits_dict = {key: edits_dict[key] for key in ("A45G",)}

#for homology, MOI, and temperature parameters use central_edits_dict (pSBK.164) and dRT_edits_dict (pSBK.160)

#for edit placement, add keys individually 

#outside fuzzy search for internal L-tail editor on lambda + dRT fuzzy searches for L-tail stop codon editor on lambda 
L_fuzzy_search_out = "TATGACCAGCCAACGTCCGA"
R_fuzzy_search_out = "ACTTTCGCAGTAAATCCCAT"
dRT_L_fuzzy_search_out = "TATGACCAGCCAACGTCCGA" # same internal lambda editor, since used same forward primer CF349 
dRT_R_fuzzy_search_out = "ATCCCATGACACAGACAGAA" # this fuzzy search is pretty far away, so may be low quality on this end... 


#set up fastq_dict = {"fastq_name": [L_fuzzy_search_out, R_fuzzy_search_out]}
fastq_dict = {"msSBK-32-14_S14_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-15_S15_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-16_S16_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-17_S17_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-18_S18_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-19_S19_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-25_S25_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-26_S26_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-27_S27_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-28_S28_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-29_S29_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-30_S30_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-36_S36_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-37_S37_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-38_S38_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-39_S39_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-40_S40_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-41_S41_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-47_S47_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-48_S48_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-49_S49_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-50_S50_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict],
             "msSBK-32-51_S51_L001_R1_001.fastq": ["TATGACCAGCCAACGTCCGA", "ACTTTCGCAGTAAATCCCAT", central_edits_dict]}

# #IF ANALYZING CI MUT EDITING 
# #outside fuzzy search for cI mut editor on lambda 
# L_cI_out = "GGCAGTCAGGCGTTGGTGCT"
# R_cI_out = "ATCGCCAGAGAAATCTACGA"

# edits_dict_cI = {"T23A": {"L_fuzzy_inside": "ATCAATGCAT", "R_fuzzy_inside": "AAATGCTTAT"},
#                  "A52T": {"L_fuzzy_inside": "ATTGCTTGCA", "R_fuzzy_inside": "AAATTCTCAA"}}


# direct it to where the reads are
data_loc = "./data/homology_unzipped/"  

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


