from Bio import SeqIO
import fuzzysearch

def determine_sequence_edits(read, L_fuzzy_search_out, R_fuzzy_search_out, edits_dict, fuzziness=1):
    """
    This function determines whether the reads map to a 
    region of interest out of the donor, if they make to a certain region of interest in the donor, 
    and if there is an edit at a certain base position within that region.
    It can look for multiple edits at once. 

    Inputs:
    read:                       a single read (as a Seq object)
    L/R_fuzzy_search_out:       a search string outside of the edit
                                (current use: 18 bases just outside the donor on phage genome)
                                L/R_fuzzy_search_in are 10 bases in the donor
                                (the first and last edit have fuzzies that extend 
                                a few bases off the donor, since they are too close to the end)
                                surrounding each edit
    fuzziness:                  how many differences from the WT sequence are allowed.
    edits_dict:                 a dictionary where the key is the edit named (ex. A4T) and the value is
                                an internally nested dictionary with two key:value pairs: a left flanking
                                search L_fuzzy_inside (right next to the edit) and a right flanking search
                                R_fuzzy_inside (just to the right of the edit)
                                the edited nucleotide is BETWEEN these two search strings
  

    Outputs:
    does_it_map:                   boolean, true if the L_fuzzy_search_out and R_fuzzy_search_out each
                                   have 1 match (add up to 2)
    does_it_map_around_edit_dict:  a dictionary with the keys as the named edit (ex. A4T) and the value
                                   as a boolean, true if the L_fuzzy_search_in and R_fuzzy_search_in each have
                                   1 match (add up to 2)
    nt_at_edit_pos_dict:           a dictionary with the keys as the named edit (ex. A4T) and the value
                                   as the nucleotide at that position
                                   if there is no mapping, the values will be "N"

    This function processes a single read. wrapper_sequence_edits.py cycles through multiple fastqs and all reads
    in a fastq to process full sets of data.

    # TO DOS:
    - change does_it_map_around_edit to a dictionary also
    """

    L_fuzzy_match_out = fuzzysearch.find_near_matches(L_fuzzy_search_out, read, max_l_dist=fuzziness)
    R_fuzzy_match_out = fuzzysearch.find_near_matches(R_fuzzy_search_out, read, max_l_dist=fuzziness)

    # construct nt_at_edit_pos_dict as if nothing maps:
    nt_at_edit_pos_dict = {}
    does_it_map_around_edit_dict = {}
    for edit in edits_dict.keys():
        nt_at_edit_pos_dict[edit] = "N"
        does_it_map_around_edit_dict[edit] = False


    does_it_map = ((len(L_fuzzy_match_out) == 1) & (len(R_fuzzy_match_out) == 1)) 

    if does_it_map == False:
        # construct nt_at_edit_pos_dict
        return {"does_it_map": does_it_map,
                "does_it_map_around_edit": does_it_map_around_edit_dict,
                "nt_at_edit_pos_dict": nt_at_edit_pos_dict}
    else:
        #zoom in on region of interest 
        zoom_read = read[L_fuzzy_match_out[0].end:R_fuzzy_match_out[0].start] 

        ##maybe this next part goes in the wrapper? 
        #values populated by "edits" dictionary 
        for edit in edits_dict.keys():
            L_fuzzy_inside = edits_dict[edit]["L_fuzzy_inside"] 
            R_fuzzy_inside = edits_dict[edit]["R_fuzzy_inside"] 

            L_fuzzy_match_in = fuzzysearch.find_near_matches(L_fuzzy_inside, zoom_read, max_l_dist=fuzziness)
            R_fuzzy_match_in = fuzzysearch.find_near_matches(R_fuzzy_inside, zoom_read, max_l_dist=fuzziness)

            does_it_map_around_edit = ((len(L_fuzzy_match_in) == 1) & (len(R_fuzzy_match_in) == 1))

            if does_it_map_around_edit == True: 
                edit_nucleotide = zoom_read[L_fuzzy_match_in[0].end:R_fuzzy_match_in[0].start]
                does_it_map_around_edit_dict[edit] = True
                nt_at_edit_pos_dict[edit] = edit_nucleotide
        return {"does_it_map": does_it_map,
                "does_it_map_around_edit": does_it_map_around_edit_dict,
                "nt_at_edit_pos_dict": nt_at_edit_pos_dict}


                     







