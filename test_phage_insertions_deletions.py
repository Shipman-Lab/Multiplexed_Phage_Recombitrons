from phage_insertions_deletions import *


def test_amplicon_insertion_deletion():
    amplicon_wt = "CCTTTGGTAAAGGTTCTAAGCTCAGGTGAGAACATCCCTGCCTGAACATGAGAAAAAACAGGGTACTCATACTCACTTCTAAGTGACGGCTGCATACTAACCGCTTCATACATCTCGTAGATTTCTCTGGCGATTGAAGGGCTAAATTCTTCAACGCTAACTTTGAGAATTTTTGCAAG"
    amplicon_4bp_insertion = "CTTTGGTAAAGGTTCTAAGCTCAGGTGAGAACATCCCTGCCTGAACATGAGAAAAAACAGGGTACTCATACTCACTTCTAAGTGACGGCTGCATACTAACCGCTTGAAGCATACATCTCGTAGATTTCTCTGGCGATTGAAGGGCTAAATTCTTCAACGCTAACTTTGAGAATTTT"
    amplicon_34bp_insertion = "CCTTTGGTAAAGGTTCTAAGCTCAGGTGAGAACATCCCTGCCTGAACATGAGAAAAAACAGGGTACTCATACTCACTTCTAAGTGACGGCTGCATACTAACCGCTTGAAGTTCCTATTCTCTAGAAAGTATAGGAACTTCCATACATCTCGTAGATTTCTCTGGCGATTGAAGGGCTAAATTCTTCAACGCTAACTTTGAGAATTTTTGCAAG"
    amplicon_4bp_deletion = "TGGTAAAGGTTCTAAGCTCAGGTGAGAACATCCCTGCCTGAACATGAGAAAAAACAGGGTACTCATACTCACTTCTAAGTGACGGCTGCATACTAACCCATACATCTCGTAGATTTCTCTGGCGATTGAAGGGCTAAATTCTTCAACGCTAACTTTGAGAATTTT"
    amplicon_32bp_deletion = "TGGTAAAGGTTCTAAGCTCAGGTGAGAACATCCCTGCCTGAACATGAGAAAAAACAGGGTACTCATACTCCATACATCTCGTAGATTTCTCTGGCGATTGAAGGGCTAAATTCTTCAACGCTAACTTTGAGAATTTT"
    off_target = "GAGGGGAGTGAAAATTCCCCTAATTCGATGAAGATTCTTGCTCAATTGTTATCAGCTATGCGCCGACCAGAACACCTTGCCGATCAGCCAAACGTCTCTTCAGGCCACTGACTAGCGATAACTTTCCCCACAACGGAACAACTCTCATTGCATGGGATCATTGGGTACTGTGGGT"

    # left & right flank should be out of any deletion regions
    left_flank = "AAAACAGGGTACTCATACTC"
    right_flank = "CATACATCTCGTAGATTTCT"

    insert_4bp = "GAAG"
    insert_34bp = "GAAGTTCCTATTCTCTAGAAAGTATAGGAACTTC"
    delete_4bp = "GCTT"
    delete_32bp = "ACTTCTAAGTGACGGCTGCATACTAACCGCTT"

    assert {"does_it_map": True, "contain_middle_search": False} == amplicon_insertion_deletion(amplicon_wt,
                                                                                     left_flank,
                                                                                     right_flank,
                                                                                     insert_34bp,
                                                                                     fuzziness=4)

    assert {"does_it_map": True, "contain_middle_search": True} == amplicon_insertion_deletion(amplicon_34bp_insertion,
                                                                                     left_flank,
                                                                                     right_flank,
                                                                                     insert_34bp,
                                                                                     fuzziness=4)


    assert {"does_it_map": True, "contain_middle_search": True} == amplicon_insertion_deletion(amplicon_4bp_insertion,
                                                                                     left_flank,
                                                                                     right_flank,
                                                                                     insert_4bp,
                                                                                     fuzziness=0)

    assert {"does_it_map": True, "contain_middle_search": False} == amplicon_insertion_deletion(amplicon_4bp_deletion,
                                                                                     left_flank,
                                                                                     right_flank,
                                                                                     delete_4bp,
                                                                                     fuzziness=0)
    
    assert {"does_it_map": True, "contain_middle_search": False} == amplicon_insertion_deletion(amplicon_32bp_deletion,
                                                                                     left_flank,
                                                                                     right_flank,
                                                                                     delete_32bp,
                                                                                     fuzziness=4)

    assert {"does_it_map": False, "contain_middle_search": False} == amplicon_insertion_deletion(off_target,
                                                                                     left_flank,
                                                                                     right_flank,
                                                                                     insert_34bp,
                                                                                     fuzziness=4)