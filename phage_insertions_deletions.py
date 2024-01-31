import fuzzysearch
import numpy as np

def amplicon_insertion_deletion(read, left_flank, right_flank, middle_search_seq, fuzziness):
    """
    Looks for insertions and deletions in an amplicon sequencing read

    Inputs:
    read:                 sequence as a string
    left_flank:           left search string (20 bp)
    right_flank:          right search string (20 bp)

    Outputs:
    does_it_map:              bool for if the read maps to the right region
    contain_middle_search:    bool for if the read contains the insertion or wt site
                              True (if insertion case): insertion occurred
                              False (if deletion case): deletion occurred
    """
    [mapped, L_fuzzy_end_index, R_fuzzy_start_index] = does_it_map(read, left_flank, right_flank)
    if mapped == False:
        return {"does_it_map": False,
                "contain_middle_search": np.NaN}
    if L_fuzzy_end_index > R_fuzzy_start_index:
        return {"does_it_map": False,
                "contain_middle_search": np.NaN}

    chopped_read = chop_read(read, L_fuzzy_end_index, R_fuzzy_start_index)
    if len(middle_search_seq) > 20:
        search_half = int(len(middle_search_seq)/2)
        middle_search_seq = middle_search_seq[search_half-10:search_half+10]

    # import pdb
    # pdb.set_trace()
    match = fuzzysearch.find_near_matches(middle_search_seq, chopped_read, max_l_dist=int(fuzziness))

    if len(match) == 1:
        return {"does_it_map": mapped,
                "contain_middle_search": True}
    elif len(match) == 0:
        return {"does_it_map": mapped,
                "contain_middle_search": False}
    else:
        return {"does_it_map": mapped,
                "contain_middle_search": "Error"}

def shotgun_insertion_deletion(read, left_flank, right_flank, middle_search_seq, fuzziness):
    """
    Looks for insertions and deletions in an shotgun sequencing read

    Inputs:
    read:                 sequence as a string
    left_flank:           left search string (20 bp)
    right_flank:          right search string (20 bp)

    Outputs:
    does_it_map:              bool for if the read maps to the right region
    contain_middle_search:    bool for if the read contains the insertion or wt site
                              True (if insertion case): insertion occurred
                              False (if deletion case): deletion occurred
    """

    L_match = fuzzysearch.find_near_matches(left_flank, read, max_l_dist=5)
    if (len(L_match) == 1):
        import pdb
        pdb.set_trace()
        read = read[L_match[0].end:]

    R_match = fuzzysearch.find_near_matches(right_flank, read, max_l_dist=5)
    if (len(R_match) == 1):
        import pdb
        pdb.set_trace()
        read = read[:R_match[0].start]

    len_required_for_search = 20
    if ((len(R_match) == 1) | (len(L_match) == 1)) & len(read) >= len_required_for_search:
        # look for inserted sequence if both flanking regions are there:
        if (len(R_match) == 1) & (len(L_match) == 1):
            middle_match = fuzzysearch.find_near_matches(middle_search_seq, read, max_l_dist=fuzziness)
        # look for inserted sequence if only left matches
        elif (len(R_match) == 0) & (len(L_match) == 1):
            middle_match = fuzzysearch.find_near_matches(middle_search_seq[0:len_required_for_search],
                                                         read, max_l_dist=fuzziness)
        # look for inserted sequence if only right matches
        elif (len(R_match) == 1) & (len(L_match) == 0):
            middle_match = fuzzysearch.find_near_matches(middle_search_seq[len_required_for_search-1:],
                                                         read, max_l_dist=fuzziness)
        if len(middle_match) == 1:
            return {"does_it_map": True,
                    "contain_middle_search": True}
        else:
            return {"does_it_map": True,
                   "contain_middle_search": False}
    else:
        return {"does_it_map": False,
                "contain_middle_search": False}


def does_it_map(read, left_flank, right_flank):
    """
    This function takes a single read & searches for sequences
    near, but not on the barcoded region using a quite fuzzy search

    Inputs:
    read:                 sequence as a string
    L_fuzzy_string:       left search string (20 bp)
    R_fuzzy_string:       right search string (20 bp)

    Outputs:
    does_it_map:          bool for if the read maps to the right region
    L_fuzzy_end_index:    index at which the left fuzzy string match ends
                          (for cutting up the read later) 
    R_fuzzy_start_index:  index at which the right fuzzy string match starts
                          (for cutting up the read later)  
    """
    L_match = fuzzysearch.find_near_matches(left_flank, read, max_l_dist=5)
    R_match = fuzzysearch.find_near_matches(right_flank, read, max_l_dist=5)

    if (len(L_match) != 1) | (len(R_match) != 1):
        does_it_map = False
        L_fuzzy_end_index = np.NaN
        R_fuzzy_start_index = np.NaN
    else:
        does_it_map = True
        L_fuzzy_end_index = L_match[0].end - 1 # weird artifact of fuzzysearch puts the match 1 up
        R_fuzzy_start_index = R_match[0].start

    return [does_it_map, L_fuzzy_end_index, R_fuzzy_start_index]


def chop_read(read, L_flank_end_index, R_flank_start_index):
    """
    Takes a read and chops it based on the L & R fuzzy matches
    """
    chopped_read = read[L_flank_end_index+1:R_flank_start_index]
    return chopped_read