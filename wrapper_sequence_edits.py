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
stop_edits_dict = {key: edits_dict[key] for key in ("A45G",)}

#for homology parameters, use central_edits_dict (pSBK.164) and stop_edits_dict (pSBK.160)
#for MOI, and temperature parameters, use stop_edits_dict (pSBK.160, pSBK.148)

#for edit placement, add keys individually 

#outside fuzzy search for internal L-tail editor on lambda + fuzzy searches for L-tail stop codon editor on lambda 
L_fuzzy_search_out = "TATGACCAGCCAACGTCCGA" #this works for edit placement, recode, and homology (forward CF349, reverse lambda_L_genoMS_R)
R_fuzzy_search_out = "ACTTTCGCAGTAAATCCCAT" #this works for edit placement, recode, and homology (except for dRT control)
stop_L_fuzzy_search_out = "TATGACCAGCCAACGTCCGA" # this works for stop codon edit pSBK.160 (dRT control) (forward CF349, reverse lambda_L_genoMS_R) 
stop_R_fuzzy_search_out = "ATCCCATGACACAGACAGAA" # this fuzzy search is pretty far away from the start, so may be low quality on this end... 
MOI_temp_L_fuzzy_search_out = "TGTAAGTTCCGCAATAACGT" # this works for stop codon edit pSBK.148 (lambda_L_genoMS_F, lambda_L_genoMS_R)
MOI_temp_R_fuzzy_search_out = "CGATGTGCGCCAGCGGAGTC" # this works for stop codon edit pSBK.148

#set up fastq_dict = {"fastq_name": [L_fuzzy_search_out, R_fuzzy_search_out]}
# fastq_dict = {"msSBK-32-63_S63_L001_R1_001.fastq":[MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
#              "msSBK-32-64_S64_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
#              "msSBK-32-65_S65_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
#              "msSBK-32-66_S66_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
fastq_dict = {"msSBK-33-89_S88_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-90_S89_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-91_S90_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-92_S91_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-93_S92_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-94_S93_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-95_S94_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-96_S95_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-97_S96_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-98_S97_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-99_S98_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-100_S99_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-101_S100_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-102_S101_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-103_S102_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-104_S103_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-105_S104_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-106_S105_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-107_S106_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-108_S107_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-109_S108_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-110_S109_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-111_S110_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-112_S111_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-113_S112_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-114_S113_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-115_S114_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-116_S115_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-117_S116_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-118_S117_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-119_S118_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-120_S119_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-121_S120_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-122_S121_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-123_S122_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],
             "msSBK-33-124_S123_L001_R1_001.fastq": [MOI_temp_L_fuzzy_search_out, MOI_temp_R_fuzzy_search_out, stop_edits_dict],}

             # "msSBK-32-64_S64_L001_R1_001.fastq": [stop_L_fuzzy_search_out, stop_R_fuzzy_search_out ,central_edits_dict],
             # "msSBK-32-65_S65_L001_R1_001.fastq": [stop_L_fuzzy_search_out, stop_R_fuzzy_search_out ,central_edits_dict],
             # "msSBK-32-66_S66_L001_R1_001.fastq": [stop_L_fuzzy_search_out, stop_R_fuzzy_search_out ,central_edits_dict],
             # "msSBK-32-67_S67_L001_R1_001.fastq": [stop_L_fuzzy_search_out, stop_R_fuzzy_search_out ,central_edits_dict],
             # "msSBK-32-68_S68_L001_R1_001.fastq": [stop_L_fuzzy_search_out, stop_R_fuzzy_search_out ,central_edits_dict],
             # "msSBK-32-69_S69_L001_R1_001.fastq": [stop_L_fuzzy_search_out, stop_R_fuzzy_search_out ,central_edits_dict],
             # "msSBK-32-70_S70_L001_R1_001.fastq": [stop_L_fuzzy_search_out, stop_R_fuzzy_search_out ,central_edits_dict],
             # "msSBK-32-76_S76_L001_R1_001.fastq": [stop_L_fuzzy_search_out, stop_R_fuzzy_search_out ,central_edits_dict],
             # "msSBK-32-77_S77_L001_R1_001.fastq": [stop_L_fuzzy_search_out, stop_R_fuzzy_search_out ,central_edits_dict],
             # "msSBK-32-78_S78_L001_R1_001.fastq": [stop_L_fuzzy_search_out, stop_R_fuzzy_search_out ,central_edits_dict],
             # "msSBK-32-79_S79_L001_R1_001.fastq": [stop_L_fuzzy_search_out, stop_R_fuzzy_search_out ,central_edits_dict],
             # "msSBK-32-80_S80_L001_R1_001.fastq": [stop_L_fuzzy_search_out, stop_R_fuzzy_search_out ,central_edits_dict],
             # "msSBK-32-81_S81_L001_R1_001.fastq": [stop_L_fuzzy_search_out, stop_R_fuzzy_search_out ,central_edits_dict],
             # "msSBK-32-82_S87_L001_R1_001.fastq": [stop_L_fuzzy_search_out, stop_R_fuzzy_search_out ,central_edits_dict],
             # "msSBK-32-83_S88_L001_R1_001.fastq": [stop_L_fuzzy_search_out, stop_R_fuzzy_search_out ,central_edits_dict],

# #IF ANALYZING CI MUT EDITING 
# #outside fuzzy search for cI mut editor on lambda 
# L_cI_out = "GGCAGTCAGGCGTTGGTGCT"
# R_cI_out = "ATCGCCAGAGAAATCTACGA"

# edits_dict_cI = {"T23A": {"L_fuzzy_inside": "ATCAATGCAT", "R_fuzzy_inside": "AAATGCTTAT"},
#                  "A52T": {"L_fuzzy_inside": "ATTGCTTGCA", "R_fuzzy_inside": "AAATTCTCAA"}}


# direct it to where the reads are
data_loc = "./data/MOI_temp_fastqs_unzipped/"  

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


