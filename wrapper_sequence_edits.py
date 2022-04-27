from Bio import SeqIO
from determine_sequence_edits import determine_sequence_edits
import pandas as pd

# define how to analyze each fastq 
# manual input fastq names and search variables 
#fastq_dict = {"msSBK-27-29_S29_L001_R1_001.fastq": [L_fuzzy_search, R_fuzzy_search, L_to_R_length, edit_position],
              #"fastq_name2": [L_fuzzy_search, R_fuzzy_search, L_to_R_length, edit_position]}

# (for testing)
#fastq_dict = {"msSBK-27-29_S29_L001_R1_001.fastq": 0}
             

# direct it to where the reads are
data_loc = "./data/"  

# loop through all the fastqs
# handle is the fastq 
for fastq_name in fastq_dict.keys():
    print(fastq_name)
    with open(data_loc+fastq_name) as handle:
        for record in SeqIO.parse(handle, "fastq"):
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