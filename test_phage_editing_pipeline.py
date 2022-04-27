import pytest
from determine_sequence_edits import determine_sequence_edits

def test_determine_sequence_edits():
    # testing_matrix = testing_matrix()

    WT_sequence = "CAGCCAACGTCCGATATCACGAAGGATAAATGCAGCAAATGCCTGAGCGGTTGTAAGTTCCGCAATAACGTCGGCAACTTTGGCGGCTTCCTTTCCATTAACAAACTTTCG"\
                  "CAGTAAATCCCATGACACAGACAGAATCAGCGATTCTGGCGCACGCCCGGCGATGTGCGCCAGCGGAGTCGTGCGGCTTCGTGGTAAGCACGCCGGAGGGGGAAAGATATTTC"\
                  "CCCTGCGTGAATATCTCC"

    #edit at position 44 (from pSBK.148) A to G
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


    #edit position at 2, G to A 
    edit_2_sequence = "CAGCCAACGTCCGATATCACGAAGGATAAATGCAGCAAATGCCTGAGCGGTTGTAAGTTCCGCAATAACGTCGACAACTTTGGCGGCTTCCTTTCCATTAACAAACTTTCG"\
                      "CAGTAAATCCCATGACACAGACAGAATCAGCGATTCTGGCGCACGCCCGGCGATGTGCGCCAGCGGAGTCGTGCGGCTTCGTGGTAAGCACGCCGGAGGGGGAAAGATATTTC"\
                      "CCCTGCGTGAATATCTCC"

    # 3 L_fuzzy_search mutations
    too_fuzzed_sequence = "CAGCCAACGTCCGATATCACGAAGGATAAATGCAGCAAATGCCTGAGCGGTTTTAAGTTCCGCTATACCGTCGGCAACTTTGGCGGCTTCCTTTCCATTAACAAACTTTCG"\
                          "CAGTAAATCCCATGACACAGACAGAATCAGCGATTCTGGCGCACGCCCGGCGATGTGCGCCAGCGGAGTCGTGCGGCTTCGTGGTAAGCACGCCGGAGGGGGAAAGATATTTC"\
                          "CCCTGCGTGAATATCTCC"


    wt_output = determine_sequence_edits(read=WT_sequence, L_fuzzy_search="TGTAAGTTCCGCAATAACGT", 
                                         R_fuzzy_search="CGATGTGCGCCAGCGGAGTC", L_to_R_length=90,
                                         edit_position=45, fuzziness=1)

    edit_output = determine_sequence_edits(read=edit_sequence, L_fuzzy_search="TGTAAGTTCCGCAATAACGT", 
                                           R_fuzzy_search="CGATGTGCGCCAGCGGAGTC", L_to_R_length=90,
                                           edit_position=45, fuzziness=1)

    off_target_output = determine_sequence_edits(read=off_target_sequence, L_fuzzy_search="TGTAAGTTCCGCAATAACGT", 
                                                   R_fuzzy_search="CGATGTGCGCCAGCGGAGTC", L_to_R_length=90,
                                                   edit_position=45, fuzziness=1)

    deletion_output = determine_sequence_edits(read=deletion_sequence, L_fuzzy_search="TGTAAGTTCCGCAATAACGT", 
                                                 R_fuzzy_search="CGATGTGCGCCAGCGGAGTC", L_to_R_length=90,
                                                 edit_position=45, fuzziness=1)

    edit_2_output = determine_sequence_edits(read=edit_2_sequence, L_fuzzy_search="TGTAAGTTCCGCAATAACGT", 
                                               R_fuzzy_search="CGATGTGCGCCAGCGGAGTC", L_to_R_length=90,
                                               edit_position=3, fuzziness=1)

    too_fuzzed_output = determine_sequence_edits(read=too_fuzzed_sequence, L_fuzzy_search="TGTAAGTTCCGCAATAACGT", 
                                                   R_fuzzy_search="CGATGTGCGCCAGCGGAGTC", L_to_R_length=90,
                                                   edit_position=45, fuzziness=1)



    assert wt_output == {"does_it_map": True, "correct_length": True, "nt_at_edit_pos": "A"}
    assert edit_output == {"does_it_map": True, "correct_length": True, "nt_at_edit_pos": "G"}
    assert off_target_output == {"does_it_map": False, "correct_length": False, "nt_at_edit_pos": "N"}
    assert deletion_output == {"does_it_map": True, "correct_length": False, "nt_at_edit_pos": "N"}
    assert edit_2_output == {"does_it_map": True, "correct_length": True, "nt_at_edit_pos": "A"}
    assert too_fuzzed_output == {"does_it_map": False, "correct_length": False, "nt_at_edit_pos": "N"}
