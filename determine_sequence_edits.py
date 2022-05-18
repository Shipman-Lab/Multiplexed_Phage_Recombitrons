from Bio import SeqIO
import fuzzysearch

def determine_sequence_edits(read, L_fuzzy_search, R_fuzzy_search, L_to_R_length, edit_position_array, fuzziness=1):
    """
    This function is determining whether the reads map to a 
    region of interest and if there is an edit at a certain base position within that region.
    It can be extended to look for multipe edits, but it doesn't do that currently. 

    Inputs: Read is a single read, L_fuzzy_search is the left side of read that contains the WT genome, 
    R_fuzzy_search is the right side of read that contains the WT genome, 
    fuzziness is how many differences from the WT sequence are allowed,
    L_to_R_length is how long the region of interest is, 
    edit_position_array is the base position of where to look for edits, 
            
            example -- 
            for L-tail edits with pSBK.148:
            L_fuzzy_search: TGTAAGTTCCGCAATAACGT
            R_fuzzy_search: CGATGTGCGCCAGCGGAGTC
            L_to_R_length: 90
            edit_position: 45 (defined from geneious, so edit_position-1)
  

    Outputs: does_it_map is true if the L_fuzzy_search and R_fuzzy_search each have 1 match (add up to 2)
             correct_length is true if the L_to_R_length fits = the search length 
             nt_at_edit_pos is the nucleotide at the edit position

             use wrapper script to cycle through all reads and count edits 

    """

    L_fuzzy_match = fuzzysearch.find_near_matches(L_fuzzy_search, read, max_l_dist=fuzziness)
    R_fuzzy_match = fuzzysearch.find_near_matches(R_fuzzy_search, read, max_l_dist=fuzziness)

    # import pdb
    # pdb.set_trace()

    does_it_map = ((len(L_fuzzy_match) == 1) & (len(R_fuzzy_match) == 1)) 

    if does_it_map == False:
        nt_at_edit_pos = "N"
        correct_length = False
        return {"does_it_map": does_it_map,
                "correct_length": correct_length,
                "nt_at_edit_pos": nt_at_edit_pos}
    else:
        correct_length = (L_to_R_length == R_fuzzy_match[0].start - L_fuzzy_match[0].end) 
        if correct_length == True: 
            chopped_read = read[L_fuzzy_match[0].end:R_fuzzy_match[0].start]
            nt_array = []
            for base_position in edit_position_array: 
                nt_at_edit_pos = chopped_read[base_position-1]
                nt_array.append(nt_at_edit_pos)
            return {"does_it_map": does_it_map,
                    "correct_length": correct_length,
                    "nt_at_edit_pos": nt_array}
        else:
            nt_at_edit_pos = "N"
            correct_length = False
            return {"does_it_map": does_it_map,
                    "correct_length": correct_length,
                    "nt_at_edit_pos": nt_at_edit_pos}


                    







