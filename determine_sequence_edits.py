from Bio import SeqIO
import fuzzysearch

def determine_sequence_edits(read, L_fuzzy_search_out, R_fuzzy_search_out, L_fuzzy_search_in, R_fuzzy_search_in, fuzziness=1):
    """
    This function determines whether the reads map to a 
    region of interest out of the donor, if they make to a certain region of interest in the donor, 
    and if there is an edit at a certain base position within that region.
    It can look for multiple edits at once. 

    Inputs: Read is a single read, L/R_fuzzy_search_out 18 bases just outside the donor, 
    L/R_fuzzy_search_in are 10 bases in the donor (the first and last edit have fuzzies that extend 
    a few bases off the donor, since they are too close to the end) surrounding each edit,
    fuzziness is how many differences from the WT sequence are allowed.
  

    Outputs: does_it_map is true if the L_fuzzy_search_out and R_fuzzy_search_out each have 1 match (add up to 2)
             does_it_map_around_edit is true if the L_fuzzy_search_in and R_fuzzy_search_in each have 1 match (add up to 2)
             nt_at_edit_pos is the nucleotide at the edit position
             for multiple edits, this nt_at_edit_pos is a dictionary -- nt_dict that has keys for each edit and values for the base at the edit site 

             use wrapper script to cycle through all reads and count edits 

    """

    L_fuzzy_match_out = fuzzysearch.find_near_matches(L_fuzzy_search_out, read, max_l_dist=fuzziness)
    R_fuzzy_match_out = fuzzysearch.find_near_matches(R_fuzzy_search_out, read, max_l_dist=fuzziness)

    # import pdb
    # pdb.set_trace()

    does_it_map = ((len(L_fuzzy_match_out) == 1) & (len(R_fuzzy_match_out) == 1)) 

    if does_it_map == False:
        nt_at_edit_pos = "N"
        does_it_map_around_edit = False
        return {"does_it_map": does_it_map,
                "does_it_map_around_edit": does_it_map_around_edit,
                "nt_at_edit_pos": nt_at_edit_pos}
    else:
        #zoom in on region of interest 
        zoom_read = read[L_fuzzy_match_out[0].end:R_fuzzy_match_out[0].start] 

        ##maybe this next part goes in the wrapper? 
        #values populated by "edits" dictionary 
        
        nt_dict = {}

        for edit in edits_dict.keys():
            L_fuzzy_inside = edits_dict[edit]["L_fuzzy_inside"] 
            R_fuzzy_inside = edits_dict[edit]["R_fuzzy_inside"] 

            L_fuzzy_match_in = fuzzysearch.find_near_matches(L_fuzzy_inside, zoom_read, max_l_dist=fuzziness)
            R_fuzzy_match_in = fuzzysearch.find_near_matches(R_fuzzy_inside, zoom_read, max_l_dist=fuzziness)


            does_it_map_around_edit = ((len(L_fuzzy_match_in) == 1) & (len(R_fuzzy_match_in) == 1))

            if does_it_map_around_edit == True:  
                chopped_read = read[L_fuzzy_match_in[0].end:R_fuzzy_match_in[0].start]
            
                for base_position in edit_position_array: 
                    nt_at_edit_pos = chopped_read[base_position-1]
                    nt_dict[edit] = nt_at_edit_pos
                return {"does_it_map": does_it_map,
                        "does_it_map_around_edit": does_it_map_around_edit,
                        "nt_at_edit_pos": nt_dict}
            else:
                nt_at_edit_pos = "N"
                does_it_map_around_edit = False
                return {"does_it_map": does_it_map,
                        "does_it_map_around_edit": does_it_map_around_edit,
                        "nt_at_edit_pos": nt_at_edit_pos}


                     







