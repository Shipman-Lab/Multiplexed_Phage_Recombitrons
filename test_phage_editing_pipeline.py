import pytest
from determine_sequence_edits import determine_sequence_edits

"""
NOTES TO SELF:
- add a test case with multiple edits where one edit maps and the other edit doesn't
- add a test case with single edit where there is a deletion between the inner search
  or insertion
"""

def test_determine_sequence_edits():
    # testing_matrix = testing_matrix()

    WT_sequence = "CAGCCAACGTCCGATATCACGAAGGATAAATGCAGCAAATGCCTGAGCGGTTGTAAGTTCCGCAATAACGTCGGCAACTTTGGCGGCTTCCTTTCCATTAACAAACTTTCG"\
                  "CAGTAAATCCCATGACACAGACAGAATCAGCGATTCTGGCGCACGCCCGGCGATGTGCGCCAGCGGAGTCGTGCGGCTTCGTGGTAAGCACGCCGGAGGGGGAAAGATATTTC"\
                  "CCCTGCGTGAATATCTCC"

    #edit at position 44 (from pSBK.148) A to G (0-44)
    edit_sequence = "CAGCCAACGTCCGATATCACGAAGGATAAATGCAGCAAATGCCTGAGCGGTTGTAAGTTCCGCAATAACGTCGGCAACTTTGGCGGCTTCCTTTCCATTAACAAACTTTCGCAGTGA"\
                    "ATCCCATGACACAGACAGAATCAGCGATTCTGGCGCACGCCCGGCGATGTGCGCCAGCGGAGTCGTGCGGCTTCGTGGTAAGCACGCCGGAGGGGGAAAGATATTTCCCCTGCGTGAATATCTCC"

    #random L-tail sequence, same length 
    off_target_sequence = "CTGAATGGCAAAGGCACCAGTACGCGCCCCACGCTGACGGTTTCTAACCTGTACGGTATGGTCACCGGGATGGCGGAAGATATGCAGAGTCTGGTCGGCGG"\
                          "AACGGTGGTCCGGCGTAAGGTTTACGCCCGTTTTCTGGATGCGGTGAACTTCGTCAACGGAAACAGTTACGCCGATCCGGAGCAGGAGGTGATCAGCCGC"\
                          "TGGCGCATTGAGCAGTGCAGCGAACTGAGCGCGGTGAGTGC"

    #large deletion in precise region 
    deletion_sequence = "CAGCCAACGTCCGATATCACGAAGGATAAATGCAGCAAATGCCTGAGCGGTTGTAAGTTCCGCAATAACGTCGGCAACTTTGGCGGCT"\
                        "TCCTTTCCATTAACAAACTTTCGGCCCGGCGATGTGCGCCAGCGGAGTCGTGCGGCTTCGTGGTAAGCACGCCGGAGGGGGAAAGATATTTC"\
                        "CCCTGCGTGAATATCTCC"


    #edit position at 2 (0,1,2), G to A 
    edit_2_sequence = "CAGCCAACGTCCGATATCACGAAGGATAAATGCAGCAAATGCCTGAGCGGTTGTAAGTTCCGCAATAACGTCGACAACTTTGGCGGCTTCCTTTCCATTAACAAACTTTCG"\
                      "CAGTAAATCCCATGACACAGACAGAATCAGCGATTCTGGCGCACGCCCGGCGATGTGCGCCAGCGGAGTCGTGCGGCTTCGTGGTAAGCACGCCGGAGGGGGAAAGATATTTC"\
                      "CCCTGCGTGAATATCTCC"

    # 3 L_fuzzy_search mutations
    too_fuzzed_sequence = "CAGCCAACGTCCGATATCACGAAGGATAAATGCAGCAAATGCCTGAGCGGTTTTAAGTTCCGCTATACCGTCGGCAACTTTGGCGGCTTCCTTTCCATTAACAAACTTTCG"\
                          "CAGTAAATCCCATGACACAGACAGAATCAGCGATTCTGGCGCACGCCCGGCGATGTGCGCCAGCGGAGTCGTGCGGCTTCGTGGTAAGCACGCCGGAGGGGGAAAGATATTTC"\
                          "CCCTGCGTGAATATCTCC"

    # 2 edits GG to TT at position 1, 2 (0,1,2)
    two_edits_sequence = "CAGCCAACGTCCGATATCACGAAGGATAAATGCAGCAAATGCCTGAGCGGTTGTAAGTTCCGCAATAACGTCTTCAACTTTGGCGGCTTCCTTTCCATTAACAAACTTTCG"\
                  "CAGTAAATCCCATGACACAGACAGAATCAGCGATTCTGGCGCACGCCCGGCGATGTGCGCCAGCGGAGTCGTGCGGCTTCGTGGTAAGCACGCCGGAGGGGGAAAGATATTTC"\
                  "CCCTGCGTGAATATCTCC"

    dict1 = {"A44G": {"L_fuzzy_inside": "CTTTCGCAGT", "R_fuzzy_inside": "AATCCCATGA"}}
    dict2 = {"G2A": {"L_fuzzy_inside": "CAATAACGTCG", "R_fuzzy_inside": "CAACTTTGG"}}
    dict3 = {"G1T": {"L_fuzzy_inside": "CAATAACGTC", "R_fuzzy_inside": "NCAACTTTGG"},
             "G2T": {"L_fuzzy_inside": "AATAACGTCN", "R_fuzzy_inside": "CAACTTTGGC"}} 

    wt_output = determine_sequence_edits(read=WT_sequence, L_fuzzy_search_out="TGTAAGTTCCGCAATAACGT", 
                                         R_fuzzy_search_out="CGATGTGCGCCAGCGGAGTC", edits_dict=dict1,
                                         fuzziness=1)

    edit_output = determine_sequence_edits(read=edit_sequence, L_fuzzy_search_out="TGTAAGTTCCGCAATAACGT", 
                                           R_fuzzy_search_out="CGATGTGCGCCAGCGGAGTC", edits_dict=dict1,
                                           fuzziness=1)

    off_target_output = determine_sequence_edits(read=off_target_sequence, L_fuzzy_search_out="TGTAAGTTCCGCAATAACGT", 
                                                 R_fuzzy_search_out="CGATGTGCGCCAGCGGAGTC", edits_dict=dict1,
                                                 fuzziness=1)

    deletion_output = determine_sequence_edits(read=deletion_sequence, L_fuzzy_search_out="TGTAAGTTCCGCAATAACGT", 
                                               R_fuzzy_search_out="CGATGTGCGCCAGCGGAGTC", edits_dict=dict1,
                                               fuzziness=1)

    # import pdb
    # pdb.set_trace()
    edit_2_output = determine_sequence_edits(read=edit_2_sequence, L_fuzzy_search_out="GAGCGGTTGTAAGTTCCG", 
                                             R_fuzzy_search_out="CGATGTGCGCCAGCGGAGTC", edits_dict=dict2,
                                             fuzziness=1)

    # too_fuzzed_output = determine_sequence_edits(read=too_fuzzed_sequence, L_fuzzy_search_out="TGTAAGTTCCGCAATAACGT", 
    #                                              R_fuzzy_search_out="CGATGTGCGCCAGCGGAGTC", edits_dict=dict1,
    #                                              fuzziness=1)

    two_edits_output = determine_sequence_edits(read=two_edits_sequence, L_fuzzy_search_out="CCTGAGCGGTTGTAAGTTCC", 
                                                R_fuzzy_search_out="CGATGTGCGCCAGCGGAGTC", edits_dict=dict3,
                                                fuzziness=1)
    assert two_edits_output == {"does_it_map": True, "does_it_map_around_edit": {"G1T": True, "G2T": True},
                                "nt_at_edit_pos_dict": {"G1T": "T", "G2T": "T"}}
    assert wt_output == {"does_it_map": True, "does_it_map_around_edit": {"A44G": True}, "nt_at_edit_pos_dict": {"A44G": "A"}}
    assert edit_output == {"does_it_map": True, "does_it_map_around_edit": {"A44G": True}, "nt_at_edit_pos_dict": {"A44G": "G"}}
    assert off_target_output == {"does_it_map": False, "does_it_map_around_edit": {"A44G": False}, "nt_at_edit_pos_dict": {"A44G": "N"}}
    assert deletion_output == {"does_it_map": True, "does_it_map_around_edit": {"A44G": False}, "nt_at_edit_pos_dict": {"A44G": "N"}}
    assert edit_2_output == {"does_it_map": True, "does_it_map_around_edit": {"G2A": True}, "nt_at_edit_pos_dict": {"G2A": "A"}}
    # assert too_fuzzed_output == {"does_it_map": False, "does_it_map_around_edit": True, "nt_at_edit_pos_dict": "N"}

