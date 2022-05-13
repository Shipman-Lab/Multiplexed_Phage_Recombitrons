from Bio import SeqIO
from determine_sequence_edits import determine_sequence_edits
import pandas as pd

# define how to analyze each fastq 
# manual input fastq names and search variables 
#fastq_dict = {"msSBK-27-29_S29_L001_R1_001.fastq": [L_fuzzy_search, R_fuzzy_search, L_to_R_length, edit_position],
              #"fastq_name2": [L_fuzzy_search, R_fuzzy_search, L_to_R_length, edit_position]}

# (for testing)
fastq_dict = {"msSBK-29-37_S37_L001_R1_001.fastq": ["TGTAAGTTCCGCAATAACGT", "CGATGTGCGCCAGCGGAGTC", 90, 45]}
             

# direct it to where the reads are
data_loc = "./data/"  

# loop through all the fastqs
# handle is the opened fastq file
# each read is a record in iterator
# iterator (big list that you can loop through) is made by SeqIO.parse 

for fastq_name in fastq_dict.keys():
    print(fastq_name)
    with open(data_loc+fastq_name) as handle:
        df = pd.DataFrame(index = SeqIO.index(data_loc + fastq_name, 'fastq'), columns = ["does_it_map", "correct_length", "nt_at_edit_pos"])
        for read in SeqIO.parse(handle, "fastq"):
            L_fuzzy_search = fastq_dict[fastq_name][0]
            R_fuzzy_search = fastq_dict[fastq_name][1]
            L_to_R_length = fastq_dict[fastq_name][2]
            edit_position = fastq_dict[fastq_name][3]
            output = determine_sequence_edits(str(read.seq), L_fuzzy_search, R_fuzzy_search, L_to_R_length, edit_position)
            df.loc[read.id, "does_it_map"] = output['does_it_map']
            df.loc[read.id, "correct_length"] = output['correct_length']
            df.loc[read.id, "nt_at_edit_pos"] = output['nt_at_edit_pos']
        df.to_csv(data_loc + fastq_name[:-6] + "_processed_data.csv")

            #import pdb 
            #pdb.set_trace()

#use record.seq to get sequence string
#from there, look at certain position 


#open fastq
#give matrix values 
#identify individual reads
#send each read through determine_sequence_edits.py
#put into csv: nucleotide at each edit position for each read
#loop back
#save csv 