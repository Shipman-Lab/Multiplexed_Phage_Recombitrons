import fuzzysearch

def extract_and_match(read, index, rep, L_outside, R_outside, L_inside, R_inside, wt_nt, edited_nt):
    left_outside = fuzzysearch.find_near_matches(L_outside, read, max_l_dist=2)
    right_outside = fuzzysearch.find_near_matches(R_outside, read, max_l_dist=2)
    if len(left_outside) == 1 and len(right_outside) == 1:
        left_inside = fuzzysearch.find_near_matches(L_inside, read, max_l_dist=1)
        right_inside = fuzzysearch.find_near_matches(R_inside, read, max_l_dist=1)
        if len(left_inside) == 1 and len(right_inside) == 1:
            var_nt = read[left_inside[0].end:right_inside[0].start]
            if var_nt == wt_nt:
                return 'wt'
            elif var_nt == edited_nt:
                return 'edited'
            else:
                return 'unmatched_edit_nt'
        else:
            return 'unmatched_region'
    else:
        return 'unmatched_region'