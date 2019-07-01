from __future__ import print_function

from qcat import scanner, adapters
from qcat import utils
from qcat import cli
from qcat import config
# from qcat import calibration
from qcat.scanner import get_adapter_by_name
from qcat.scanner_base import find_best_adapter_template, extract_align_sequence
from qcat.scanner_epi2me import BarcodeScannerEPI2ME

barcode_spacer = "NNNNNNNNNNNNNNNNNNNNNNNN"


def test_get_placeholder():
    template = "NNNNN"
    bc = adapters.AdapterLayout.get_placeholder_pos(template)
    assert bc.start == 0
    assert bc.end == 4
    assert bc.length == 5

    template = "AAAANNNNN"
    bc = adapters.AdapterLayout.get_placeholder_pos(template)
    assert bc.start == 4
    assert bc.end == 8
    assert bc.length == 5

    template = "NNNNNAAAA"
    bc = adapters.AdapterLayout.get_placeholder_pos(template)
    assert bc.start == 0
    assert bc.end == 4
    assert bc.length == 5

    template = ""
    bc = adapters.AdapterLayout.get_placeholder_pos(template)
    assert bc.start == -1
    assert bc.end == -1
    assert bc.length == 0

    template = "AATGTACTTCGTTCAGTTACGTATTGCTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"
    bc = adapters.AdapterLayout.get_placeholder_pos(template)
    assert bc.start == -1
    assert bc.end == -1
    assert bc.length == 0

    template = "AATGTACTTCGTTCAGTTACGTATTGCTNGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"
    bc = adapters.AdapterLayout.get_placeholder_pos(template)
    assert bc.start == 28
    assert bc.end == 28
    assert bc.length == 1

    template = "AATGTACTTCGTTCAGTTACGTATTGCTNNNNNNNNNNNNNNNNNNNNNNNNGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"
    bc = adapters.AdapterLayout.get_placeholder_pos(template)
    assert bc.start == 28
    assert bc.end == 51
    assert bc.length == 24

    template = "AATGTACTTCGTTCAGTTACGTATTGCTNNNNNNNNNNNNNNNNNNNNNNNNGTTTTCGCATTTATCGTGNNNNNNNNNNAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"
    bc1 = adapters.AdapterLayout.get_placeholder_pos(template, 0)
    bc2 = adapters.AdapterLayout.get_placeholder_pos(template, 1)
    assert bc1.start == 28
    assert bc1.end == 51
    assert bc1.length == 24

    assert bc2.start == 70
    assert bc2.end == 79
    assert bc2.length == 10


def test_barcode():
    read = "AATGTACTTCGTTCAGTTACGTATTGCTTCGATTCCGTTTGTAGTCGTCTGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCAA" \
           "ATCATAAACGCGCTGGAAGTTGCGCCGTTCAATCACCATCACCTTTCCGGCAGTAAATAAACCTTCCAGATGCCGTTTATGCGGCTTCCATTCCCACCAGCCG" \
           "CTTGCACCTTTACGAGGATGCTCAAAATCGGCTGAACGTACCGGCCCCTTATCATGAATATGCTGAATTAACTGTGCAATTTCCGCCTCATGTTCCTGCATCC" \
           "AGGCGTCTTTGTATTTCCAGCCCATTTTTTCAGGTGCCAGCATGCGGTGGCGAATAAGACGAAAGTCGCTACGCGGCATAAAGCAGGCTTCATGCGCCCAGTA" \
           "TTCCATTAATTCGCCACGCGCCAGAGACTCATCCAGCCACTGGGCAGGATAATTTCCCAGACGACTGAAAAGCACCAGATATGGACTACGGGCAACAATATTG" \
           "ATGGTATCGATTTGCAGCAAGGACATGCGGGAGATCGTTGCCGGAATATCCTCCAACGACGCTCGACGGCGGGGTTTGTTTAACAGGCCTTGTGCGGCAAGGT" \
           "GAAGATTACGCGCATCAGCAAGGGAGAGGTGCGGCAGCGACATTCATGACTCCATCAATCGAACGCTGCCGCGGCGTAACTAGTTGCCAGAAGCCAGCAAGGT" \
           "TAGTTGCGTAAGCAGTTTCGCTGGTTCATCACCTGAAAGCTGTGCGTCTACAGGCAAATACCACCAATTTTCTTCTGCAAAGGCCCGGCATTTCACCGCATCT" \
           "TTTTCAGTCATTACCAGCGTTTGCCCGGCGCTTACCAACGCACTGACATCCGCATGGTTCAAAGACTGATGATCGGCCAGCGGTACACATTTTTCCGGTTGTA" \
           "CGCCACACATCTTCAGCGTGGCAAAAAAGCGCGGCGGATGCCCAATCCCCGCCATCGCCACTACATGTTCAAGCTGAGCAACGTCACAACGCGTACCGGTACG" \
           "TAAATTCACCGCCTGACCCGGCAGCAGATGCATGGGGATTTCACCGCTGCGAGGGACACCGCCGTTGACGATTACCGCATCAACCGACTTTAAGCGCCCCGCT" \
           "CGCTCAC"

    barcode_dict = BarcodeScannerEPI2ME().detect_barcode(read)
    assert barcode_dict['barcode'].name == "barcode02"


read = "AATGTACTTCGTTCAGTTACGTATTGCTTCGATTCCGTTTGTAGTCGTCTGTGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCAA" \
       "ATCATAAACGCGCTGGAAGTTGCGCCGTTCAATCACCATCACCTTTCCGGCAGTAAATAAACCTTCCAGATGCCGTTTATGCGGCTTCCATTCCCACCAGCCG" \
       "CTTGCACCTTTACGAGGATGCTCAAAATCGGCTGAACGTACCGGCCCCTTATCATGAATATGCTGAATTAACTGTGCAATTTCCGCCTCATGTTCCTGCATCC" \
       "AGGCGTCTTTGTATTTCCAGCCCATTTTTTCAGGTGCCAGCATGCGGTGGCGAATAAGACGAAAGTCGCTACGCGGCATAAAGCAGGCTTCATGCGCCCAGTA" \
       "TTCCATTAATTCGCCACGCGCCAGAGACTCATCCAGCCACTGGGCAGGATAATTTCCCAGACGACTGAAAAGCACCAGATATGGACTACGGGCAACAATATTG" \
       "ATGGTATCGATTTGCAGCAAGGACATGCGGGAGATCGTTGCCGGAATATCCTCCAACGACGCTCGACGGCGGGGTTTGTTTAACAGGCCTTGTGCGGCAAGGT" \
       "GAAGATTACGCGCATCAGCAAGGGAGAGGTGCGGCAGCGACATTCATGACTCCATCAATCGAACGCTGCCGCGGCGTAACTAGTTGCCAGAAGCCAGCAAGGT" \
       "TAGTTGCGTAAGCAGTTTCGCTGGTTCATCACCTGAAAGCTGTGCGTCTACAGGCAAATACCACCAATTTTCTTCTGCAAAGGCCCGGCATTTCACCGCATCT" \
       "TTTTCAGTCATTACCAGCGTTTGCCCGGCGCTTACCAACGCACTGACATCCGCATGGTTCAAAGACTGATGATCGGCCAGCGGTACACATTTTTCCGGTTGTA" \
       "CGCCACACATCTTCAGCGTGGCAAAAAAGCGCGGCGGATGCCCAATCCCCGCCATCGCCACTACATGTTCAAGCTGAGCAACGTCACAACGCGTACCGGTACG" \
       "TAAATTCACCGCCTGACCCGGCAGCAGATGCATGGGGATTTCACCGCTGCGAGGGACACCGCCGTTGACGATTACCGCATCAACCGACTTTAAGCGCCCCGCT" \
       "CGCTCAC"

read_bc3_exact = \
        "AATGTACTTCGTTCAGTTACGTATTGCTGAGTCTTGTGTCCCAGTTACCAGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCAA" \
        "ATCATAAACGCGCTGGAAGTTGCGCCGTTCAATCACCATCACCTTTCCGGCAGTAAATAAACCTTCCAGATGCCGTTTATGCGGCTTCCATTCCCACCAGCCG" \
        "CTTGCACCTTTACGAGGATGCTCAAAATCGGCTGAACGTACCGGCCCCTTATCATGAATATGCTGAATTAACTGTGCAATTTCCGCCTCATGTTCCTGCATCC" \
        "AGGCGTCTTTGTATTTCCAGCCCATTTTTTCAGGTGCCAGCATGCGGTGGCGAATAAGACGAAAGTCGCTACGCGGCATAAAGCAGGCTTCATGCGCCCAGTA" \
        "TTCCATTAATTCGCCACGCGCCAGAGACTCATCCAGCCACTGGGCAGGATAATTTCCCAGACGACTGAAAAGCACCAGATATGGACTACGGGCAACAATATTG" \
        "ATGGTATCGATTTGCAGCAAGGACATGCGGGAGATCGTTGCCGGAATATCCTCCAACGACGCTCGACGGCGGGGTTTGTTTAACAGGCCTTGTGCGGCAAGGT" \
        "GAAGATTACGCGCATCAGCAAGGGAGAGGTGCGGCAGCGACATTCATGACTCCATCAATCGAACGCTGCCGCGGCGTAACTAGTTGCCAGAAGCCAGCAAGGT" \
        "TAGTTGCGTAAGCAGTTTCGCTGGTTCATCACCTGAAAGCTGTGCGTCTACAGGCAAATACCACCAATTTTCTTCTGCAAAGGCCCGGCATTTCACCGCATCT" \
        "TTTTCAGTCATTACCAGCGTTTGCCCGGCGCTTACCAACGCACTGACATCCGCATGGTTCAAAGACTGATGATCGGCCAGCGGTACACATTTTTCCGGTTGTA" \
        "CGCCACACATCTTCAGCGTGGCAAAAAAGCGCGGCGGATGCCCAATCCCCGCCATCGCCACTACATGTTCAAGCTGAGCAACGTCACAACGCGTACCGGTACG" \
        "TAAATTCACCGCCTGACCCGGCAGCAGATGCATGGGGATTTCACCGCTGCGAGGGACACCGCCGTTGACGATTACCGCATCAACCGACTTTAAGCGCCCCGCT" \
        "CGCTCAC"

read_bc3 = \
    "ATGCTCAGCAATGTACTTCGTTCAGTTACGTATTGCTGAGTCTTGTGTCCCAGTTACCAGGGTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCG" \
    "CCGCTTCAAGCGCTGGAAGTTGCGCCGTTCAATCACCATCACCTTTCCGGCAGTAAATAAACCTTCCAGATGCCGTTTATGCGGCTTCCATTCCCACCAGCCG" \
    "CTTGCACCTTTACGAGGATGCTCAAAATCGGCTGAACGTACCGGCCCCTTATCATGAATATGCTGAATTAACTGTGCAATTTCCGCCTCATGTTCCTGCATCC" \
    "AGGCGTCTTTGTATTTCCAGCCCATTTTTTCAGGTGCCAGCATGCGGTGGCGAATAAGACGAAAGTCGCTACGCGGCATAAAGCAGGCTTCATGCGCCCAGTA" \
    "TTCCATTAATTCGCCACGCGCCAGAGACTCATCCAGCCACTGGGCAGGATAATTTCCCAGACGACTGAAAAGCACCAGATATGGACTACGGGCAACAATATTG" \
    "ATGGTATCGATTTGCAGCAAGGACATGCGGGAGATCGTTGCCGGAATATCCTCCAACGACGCTCGACGGCGGGGTTTGTTTAACAGGCCTTGTGCGGCAAGGT" \
    "GAAGATTACGCGCATCAGCAAGGGAGAGGTGCGGCAGCGACATTCATGACTCCATCAATCGAACGCTGCCGCGGCGTAACTAGTTGCCAGAAGCCAGCAAGGT" \
    "TAGTTGCGTAAGCAGTTTCGCTGGTTCATCACCTGAAAGCTGTGCGTCTACAGGCAAATACCACCAATTTTCTTCTGCAAAGGCCCGGCATTTCACCGCATCT" \
    "TTTTCAGTCATTACCAGCGTTTGCCCGGCGCTTACCAACGCACTGACATCCGCATGGTTCAAAGACTGATGATCGGCCAGCGGTACACATTTTTCCGGTTGTA" \
    "CGCCACACATCTTCAGCGTGGCAAAAAAGCGCGGCGGATGCCCAATCCCCGCCATCGCCACTACATGTTCAAGCTGAGCAACGTCACAACGCGTACCGGTACG" \
    "TAAATTCACCGCCTGACCCGGCAGCAGATGCATGGGGATTTCACCGCTGCGAGGGACACCGCCGTTGACGATTACCGCATCAACCGACTTTAAGCGCCCCGCT" \
    "CGCTCAC"

read_nobc = \
    "CCGCTTCAAGCGCTGGAAGTTGCGCCGTTCAATCACCATCACCTTTCCGGCAGTAAATAAACCTTCCAGATGCCGTTTATGCGGCTTCCATTCCCACCAGCCG" \
    "CTTGCACCTTTACGAGGATGCTCAAAATCGGCTGAACGTACCGGCCCCTTATCATGAATATGCTGAATTAACTGTGCAATTTCCGCCTCATGTTCCTGCATCC" \
    "AGGCGTCTTTGTATTTCCAGCCCATTTTTTCAGGTGCCAGCATGCGGTGGCGAATAAGACGAAAGTCGCTACGCGGCATAAAGCAGGCTTCATGCGCCCAGTA" \
    "TTCCATTAATTCGCCACGCGCCAGAGACTCATCCAGCCACTGGGCAGGATAATTTCCCAGACGACTGAAAAGCACCAGATATGGACTACGGGCAACAATATTG" \
    "ATGGTATCGATTTGCAGCAAGGACATGCGGGAGATCGTTGCCGGAATATCCTCCAACGACGCTCGACGGCGGGGTTTGTTTAACAGGCCTTGTGCGGCAAGGT" \
    "GAAGATTACGCGCATCAGCAAGGGAGAGGTGCGGCAGCGACATTCATGACTCCATCAATCGAACGCTGCCGCGGCGTAACTAGTTGCCAGAAGCCAGCAAGGT" \
    "TAGTTGCGTAAGCAGTTTCGCTGGTTCATCACCTGAAAGCTGTGCGTCTACAGGCAAATACCACCAATTTTCTTCTGCAAAGGCCCGGCATTTCACCGCATCT" \
    "TTTTCAGTCATTACCAGCGTTTGCCCGGCGCTTACCAACGCACTGACATCCGCATGGTTCAAAGACTGATGATCGGCCAGCGGTACACATTTTTCCGGTTGTA" \
    "CGCCACACATCTTCAGCGTGGCAAAAAAGCGCGGCGGATGCCCAATCCCCGCCATCGCCACTACATGTTCAAGCTGAGCAACGTCACAACGCGTACCGGTACG" \
    "TAAATTCACCGCCTGACCCGGCAGCAGATGCATGGGGATTTCACCGCTGCGAGGGACACCGCCGTTGACGATTACCGCATCAACCGACTTTAAGCGCCCCGCT" \
    "CGCTCAC"

real_bc03_porechop = \
    "CTCTGTACTTCGTTCAGTTACGTATTGCTGAGTCTTGTCCCCAGTTACCGGGTTTCGCATTTATCGTGAAACGCTTTCGCGTTTCGTGCGCCAACTTCACTGG" \
    "GGAATGCCGCCGATGCCGGATCAATTCTTTACCGCCAGACCTGCTTACCAGCATGGGGGCAGCCATTGGGGCCGTTAGTATGACCGGCATCCTGTTTTCTCTC" \
    "GGTGCCAGTATAGAAGATTATCAGTGATTATTAGCGCAGATGCTGGCACAAAGCCAAGAACTCCCGTATACCACAGACAACGGATAACAGTGCAGAACACCTA" \
    "TTTCCTCACTGGATAACATGGTTGCCAGAGGGCAATGTTCTACCTGTTCTGTACACGGGGAAATACCGCGTGAAGATCTGCGTGGTTTCTCAGAGAGATCAGC" \
    "AGCAGACGAAAGGACAAATTGAGTCGAAAGGTTGTGGTGATTGATTAAACTGATGCAAAATGTTTATATTAGTGACAACCTGCAGCGGGCGGTTTGTCATTAT" \
    "GGAGCGTGAGGAATGGGTAAAGTAAAGGGGCATACCCGCGCAGAAGCGAAAAGGACAACCTGAAGTCCACACCAGTTGCTGAGGAATTATGATCGATACCATC" \
    "GAAGAATTGATTGAAGGTCCGATGGATGGCTTAAAAAGCGTGCTGCTGAACAGTACGC"

simulated_bc03 = \
    "CACTGACCGCCCGCTTTGCGATTTAATGTCTTCGTTCAGTTACGTATTGCACGAGTCAGA" \
    "TGTTCCAGTTACCAGGGTTTTCTCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCC" \
    "GCTTCAAACCAACAGTAAGGCAAATCGTGTTGTACCAGCTCTAGGGCTGGATAGTGCCAG" \
    "ATTAGTCGCATTCATCTGCCGCTGGCGCTTATTTCTAAGGCGCTAAACTTTCCAAGCAGA" \
    "ATCTATGCGCCTAAGTCGCCGTAAGCCAACCACGCCGCCCGGTACTAATCGCGGGCACTT" \
    "CAATTTCTGGAGCAGCAGGCAGAAGCACACTGGCAGGATTACAGCGTCGAGCCCAAATCC" \
    "TTCAGTCAGCCGTCAAAATTGGAGGCTTACCGCCGTGCCCGAGTCGGCAATTGTAAATTC" \
    "AACATTCTCAAATGCGTCATGGCTGAGCGATGATTAGTCGCTAGAGAAGCTTGTTTTTGT" \
    "CCTGAATGATGTTTGACACTACCGAGGTGTACATTTTTACCCGAGTCGCTAATTTTTGGC" \
    "GCAAGGTGGTAAGCGCTTTGAGGAAAGCGAGGCTTGAAGCGGTAGTTCCCGTCACGTCAG" \
    "TGGTGCTCGCGGTGGACGAAACTGCTATCGACAAAGATACAGTGAAAATGCCCTGAAGGT" \
    "CTGTTCAGGGCTCAATAATGAGGGGTAGGGTAGCCTGGCTGGTTGGCGGCGGCGTGCGCG" \
    "ACCTGTTTACTTGGCAAAAAGCCGAAAGATTTTGACGTAACCACTAACGCCACGCCCGAG" \
    "CAGGTGCGCAAACGGTTCCGTAACGCTTCTGGTGCCGTAGTCTGGTCATGTAATGTTTGG" \
    "CCCGCGATTATAGAAGCTGCCTTCCCGTGGACACCACGCCAAGGAACGTCACCGACCGCA" \
    "CGACCTCCCAACGCGGGCAAATCGGCATGTTGCTGCGCGACAAACTTTCGCCTCAGAGTA" \
    "CTCGAAGAAGACGCCCAGCGCCGCGATTTCACTATCAACAGCCTGTATTACAGCGTAGCG" \
    "GATTTTACCGTCCGTGATTACGTTGGCGGACATGAAGGATCTGAAGGACGGCGTTATCCG" \
    "TCTGATTGGTAACCCGGAAACGCGCTACCGTGAAGATCCGGTACGTATGCTGCGCGCGGT" \
    "ACGTTTTGCCGCCAAATTGGGTATGCGCATCAGCCCGGAAACCGCAGAACCGATCCCTCG" \
    "CCTCGCTACCCTGCTGAACGATATGCAAGAAGCCCCACCGGCACGCCTGTTTGAAGAATC" \
    "GCTTAAACTGCTACAAGCGGGCTACGGTTACGAAACCTATAAGCTGTTGTGAATATCATG" \
    "CGTTCCAGCCGCTGTTCCCGACCATTACCCGTGACTTCACGGAAAATGGCGACAGCCCGA" \
    "TGGAGCACATCATGTTGAACAGGTGCTGAAGAGTACCGATAGAGCGTATCCATAACGATA" \
    "TGCGCTGAACCGGCGTTATCGTATGCACTGCCCACTGCTGGAGCGGACAGTGATCGCCAG" \
    "GAAAGCGGCCTGACCTATCACGACGCTTTCGCGCTGGCGATGAACGACGTGCTGGACGAA" \
    "GCCTGCGTGTCACTGGCAATCCCGAAACGTCTGACGACATTACCCGCGATATCTGCCAGT" \
    "TGCAGTTGCGTATGTCCCGTCGTCAGGGTAAACAACCGGCTTCTCGGTGAACGTCACGGC" \
    "ACCAGAGTAAGCCTGAGCGTGCCGTGCGGACTTATGACCTGGACGGGCCTGAGTCCGAGC" \
    "AGTTTAGTGTTAAAACGCTGAACTGCAGCCAGTCTGGTGAAATGGTGGGGTCAGTTCCAG" \
    "GTTTCCGCGCCACCAGACCAATAAAAGGGATGCTCAATCGAGCGGGCATGCAGAACCGTA" \
    "ATCGGATAGTCGTCCTCGTCGGACACGCAAACGCGCCCACGTCGTGAGGGTACCGCATGA" \
    "CAGTGGCGTATATTGCCGCAGCAATCTGGCCTCTCCGCTGAGCAGACATCGGGCTGCCCT" \
    "CTGAAAGCCTTAGGCGATATCCCTGAAACCACATTCTTACCTGTGTCATGTTGGTACCGC" \
    "ACCCCACCGCTGAGGCCGCAAGTCATTAACCCGACACTTAAACGCAGCCGTGGCCGCATG" \
    "AATGTGCTCTCTTGCACCTGATAAGAGCTACTTTCAATCACACAAGCGTATCACATTGCA" \
    "GCAAGGTCGCGTCCGCAAGCTGAACGCTGGGGACCACGCACGCTGGATCTCGACATCATG" \
    "CTGTTTGGTAATGAAGTGATAAATACTGAACGCCTGACCGTTCCGCACTACGATATGAAG" \
    "AATCGTGGATTTATCGCGTGGCCGCTGTTTGAAATCGCGCCGGAGTTGGTGTTTCCTGAT" \
    "GGGGAGATGTTGCGTCAAATCTTACATACAAGAGCACTTGACAAATTAAACACAATGGTA" \
    "ATTATCTATATTGCTATAATAAATATTTCTCAAACATATCTCATCTCTTTCCTTAACATT" \
    "ACGATAATAGTTTAAACATCTCTTTGCGGAATGAATATAACTTACAAAATTAAATCTCTC" \
    "CACAACCCTAAACTTGAACCACTTGTGAACTTATTTACTTAATTTTGAGTTAAGTAAATT" \
    "TAATTTGTTCAAAAGCTTCTAAAGTGCAGTATAGGCGATGTGAAAATGTTTCTATTTGGG" \
    "TAACAAATATCATTGCTCATTAATATCTAACATTTGTTTCTCCTATTTGTGTCATGAATA" \
    "GTCATTCTTCAAACAAGTAGTAATACTGGAACTGGATAAATTATGCAAGGATTTTTATGA" \
    "TCTAAAAATTAGGTTTCGCCTGTCTGGTCTTAAGTGCTGGACGATGGTTGCAGGTACTGC" \
    "TTCTGCTGATATGGACGGCGTTCAGTTAAATATCAGTTGTGGTTGTTGATAACACTCGTG" \
    "ACTCGCAGTGACGGCGGTAACAAGGATCTGATCCTGCTGCAATCGCAACCGTTGGTGACA" \
    "TCCTTGCTGCCGTCTGATTTGACACCGTTGGTCCTCTTACGCAAAGCTAAGGGTACACCA" \
    "GTTGACTGCACAGCAAAGCACCCGAATTCAGCTGGCTAAAGACCTTCTGTTCTGTTTTCT" \
    "TACTACAGAAAGGGACTTTGCACAACGGGATGCCGTCACTATTAACAACCCCTCCGCGGT" \
    "TAATATCGCTCTACACAATATTGATGGTTCCACTATCAAACAGGTTCAAATTAACAACCC" \
    "TGGCGATGTGTATACTAAAGCCCTGGATGCAACGACAAAATCTGCTGTTTATGATTTTAA" \
    "AGCGTCTTACGTTCGTGCCGTCGCAGACCAAACAGCAACTGCTGGTTAAACTTAAAAACT" \
    "AACACTGCATACACCATACTTAGTAAGTATTACAGTCAATTCGATTGAATGAATATACAG" \
    "GGAATAATAATTTCTATTTATATTATTCCCTGTTTTAATTAACTCTATCAGGGATGCTAG" \
    "ATGTTTATCAACATACCCAACAGCTTTATGCTTCGTAACCTGTATGGCTTTTAGTTCATC" \
    "GTCTATTGCGGACATTGTCCATTTCGGGTACTCGCGTAATATATAAAAGCGATCAAAAAG" \
    "TGTCAACGTACGTCTGGAAAATAAAGGGAATAACGTTGCTTGTCCAGAGTTGGTTAGATA" \
    "CTGGCGATGACAGCTGATGTCCTGACAGTATTACAGTCCATTTTACTGCTACGCCGCCAG" \
    "TATCGCGTATTGATGCCAAACGTGGGCAAACAATCAAATTAATGTACACAGCCTATACCT" \
    "CACTGCCTAAAGACAGAGAGAGCGTGTTCTGGTTTAACGTAACTGGAAGTTCCACTTAAA" \
    "CCATGCAGAGAATGTGAGATGAAGCTAATTGGTACAACTGCATTCGCTCAACGTAATTGA" \
    "AAACTTTTCGGCTGCCCGGATGGATTGAAGGGATCCGACGGAGGCCCGTTAGCCCTGAAG" \
    "TGCTCTGTCAGCTAACAATTAAGGCGACGCTCTCATATTACAAGTACTGATCCACATTTA" \
    "CTAAGTCTCTTTTTAGCAGTTGTGATTAGAAGCTAGCGGTAAACGACTATCCGATTGATG" \
    "TGAAAATGATTGCACCATTTAGTGATGATGTCATGAAAGTCAATGGCCTTAATGGCTAAG" \
    "CGAATTCTGCAAAAGTGCATTTTTACGCCATTAATGACTTTGGTGGCGCAATTGAAGGTA" \
    "GCCGGCTGCCATCGTCTGAGGGATGCGATATTTCCGCAGGAAGCATACAGCGTGACTATT" \
    "CAATATACTAAAAATTATCATACATCTGACCCGTTCGCCACGTTTTGCGCGCTGCTGTAT" \
    "TGCAATACTGCTTTCAGTGCTGAACTCGTTGAATATGACCAACCTTCCTGATGGGGCAGA" \
    "ATGCATCTAATAGTGACTCAGCCGTACAGTGAAGGTAACCCCCTATACTATTTACGACGG" \
    "AGAATCAGTGTTTTTGTAAGCAACCAATCAACCAAAGTATTACATTTGTCGCAATTGAAG" \
    "GAAAAAGACGCCCAGGCTTGGCATTCATTTAAGCAATTTATTGTTTCATATTAATTCTCC" \
    "CGATATAAATAACATAAAACCGTTCTGCTTGCCAGGGATGAACGCTCGGCAATTCCCTCA" \
    "ATTTGACGGAAATTATCCCTCAATCTTGTCGTTATGACGTTAACGATCAACGTCTGAATA" \
    "TAGTCCGACGAGACTTTCCCAGCCTGGGTAATGAAAAATTACACCCAACTATGTATAATG" \
    "CCGTTATGGGAAAACGGCATTAATGCGGCTATGTTGTCATACAACCTCAACGGATATCAT" \
    "AGTGAAACCCCTGGTCGAAAAACGAGCATTTATGCTGCTTTAACGGTGGGATGAATTTAG" \
    "GTGCATGGCGACTGCGTGCCTCGGGCAACTACAAATGGATGACCGATTCTGGCCTATTTA" \
    "TGATTTTAAGAATCGGTATGTTCAGCGTGATATCGCCGCGCTGCGTTCTCAACTCATTCT" \
    "TGGTGAGTCTATACGACGGGCGAAGTTTTGATTCCGTCAGTATCCGAGGCATTCGTTTAT" \
    "ACAGTGACAGCCTCATGTTGCCTTGATTTCTCAGCTTTGCGCCCTATCATTCATGGCGTT" \
    "CGCCAATACCAACACTAAATTAACTAATACGCAAGGTGGCTATAAGATCTAAAAGCGACG" \
    "GTGCGGCCAGGCGCTTTCGTCATTGATGATACTGAGTCCGTCAGGAACGCAGCCAGCCTA" \
    "TGTTACCATCGGAATCCGTCGCTCAAACCGGACATCGCAACCTTTCTCATCCGTTGTTCA" \
    "AATGTTACGCCCTGGCGTTGGACGTTGGGATATTAGCGGCGGTCAGGTCCCAGAATGGAT" \
    "ATTAGTAAGCTAATCTCTATTCCAAGTAAATTGCTACTACGGCCTGAAAAACTATCTGAC" \
    "GGCTTATACCGTATTCAGATAATCGATAATAACTATACCGCCGGGGTAGGTCTAGGTCTG" \
    "AATACTTCAGTTGGTGCATTTTCTTTCGACGCAACTCATTCCAATGTTCGTATCCCGGAT" \
    "GATAAAACATACCAGGGGCAAAGTTATCGTGTTTCCGGGAACAAGTTATTCGAAGAAACA" \
    "AGACTTCCTCGCCGCCTATCGCTTTTTGGAGAATTACCTTGGTCTTAACCTGATGCACTA" \
    "ACTCTAATTGATGAAGTGAAACATCCCGAACAAGATACGGAACCGAAATCCGGGCGTAAT" \
    "TACTCACGCATCTCGAATCAGGTTACGGTCAGTATTAACCAACGGTTGAAATTTGAGAAA" \
    "AAGATTCGGGCTTGTTACTTTCCGGAAGTTGCACCACTATTGGTACCTGGGGGCTTGCCG" \
    "ACAAAATCGTAGCAAGTACTCTATTGGCTACAGTAGTCATCCTCGGCAGTTACAGAGTCA" \
    "TTGTCGGATGGAGTTAAGGGATGACACTGACGATAGCGTTTATCTTAGTTTCACCATTCC" \
    "AATCGAAAAATTACTTGTGAACAAACGTACTACCGATTCAGTTTGCAGAGTATTGATACT" \
    "GCACCAGGAATAGCGGTGGCTTTAGGTATAATAACCAACTCAACGTTAGCAGCAGTGCTA" \
    "TAGCGATAACGCTCGCGTCAGTCATAGCGTGAATACTGGCTATACGCTGGCATGAATAAA" \
    "GCCAGCAAAGATTTGAGTTATGTGGGGGTTGATCCGGCTTGAGTGAACCATGGGTCGTAA" \
    "GCTGGCAGGTTCAATTTCTGCAGACATTCGAAACAGCCGCAAACTATCTCTCGGCTCGGA" \
    "CGGTGGTTTTGTATTGCATAAGCGGTGGACTGACTTTCAGTAATGATAGTTTTAGCGACT" \
    "CCGATACACTGGCGAGAGTTCAGGCTCCAGTGCTCAAGGTGCGCGAATAATTATGGCAAC" \
    "AGTACAATCGATCGATGGGGTTATGAGTGTCAGCCCCACCGCTCCTTTATCTCCTTATCA" \
    "TGGAAAACCGTATCGCGCTGAATATCAACGATCAATCGAACGTCCTGCTTGACATAAAAG" \
    "TACCAGTGCAGTAGCTGTACCGCGTCAGGGTTCAGTCGTCTTTGCTCATTTCTAAACCGT" \
    "GCAAGGGCAATCAGCCATTATGAACATCACACGAAGTGATGGTAAAAATATTCCATTTGC" \
    "TGCAGATATTTATGATGAGCAAGGCAATGTCATTGGTGATGTTGGACAGGGTGGACAAGC" \
    "ATTTGTTCGTGGTATTGAGCAGCAGGGAAATATCAGCATTAAATGGCTCGAACAAAGTAA" \
    "ACCCGTAAGTTGTCTTGCGCATTATCAACAAAGCCCAGAAGCAGAAAAATAGCACAATCT" \
    "ATTATTCTGAATGGAACCAGGTGTCAGATTCAGTAACTACAAGGGACCTCAAATGATAAA" \
    "AACATCGCCACATAAAATAGTGACACTGATGGGAATATTATTATCACCCTCAGATTTGCA" \
    "ACGGATATTAATGAAGAGTATCAGCACTGGAACACGGACTACCCGCACATCACACACTAC" \
    "TGGTATAAACGTCACGAATGATGGCAATAATAAGCTACTACGAGACTGAATCCATAAGAT" \
    "GGGTCTGGATAAGATCGCGAATAAAACGACAGAATCTCAGGCTGATTTTAACTGGTTGCC" \
    "AGTGGTGCAGCAGTGGCCTCCGTTTGGATTGATAACTGCTTGCGGAAATGCCATTACCAT" \
    "AAAGCTCACCTAAGCTTATTATACCGCAGTCTGGTTCATCTTCGACGACAAGTAATATCG" \
    "GTATGGGTTTCAAAAACGGACTACTGATGATGCCACTTTCCTTAAACCTAACAGTGCTAA" \
    "GGGAAAAGATATGCTGGAGCACAGACGAGATGCAGCCCGATAAGGGTCTTGAAATGACCG" \
    "TTGCGCAACGTGAACAGATGCAGGGCAAGGCGTACCGGGAAGCATGCTGGGCACGTTTCA" \
    "TTTTTCTATAAGTAATAAGATAATTACATGCACTGAAAAATTAGTGTTCATTAGGAGTAT" \
    "CACACAGGTCGTTTTACTTGCCTCACGTTCAGCGTACTCTCTTAGTCCAAACAGACAGCT" \
    "TGGACTGACCGTGGTGTTGAGCAATGGGTACTTGTGTACCGCACATTATTAACTACGTCG" \
    "TGACAATATTTATGTTGTTGATTTTGGGGATGTATATATTTCCGACCATACAATGCCAAT" \
    "GACCAAAGTAAAAACATTCAAACACAAATTCAAAGACTGTGCGGGTATCCCCAATAAAGG" \
    "AGACCAAATAAAATTAACTTAGCGAGCCACATGCTAGATGAACTGCTAATGACGGTGCGG" \
    "GGTTTGCAAATGGTTCCACAGCCGCAGATCAAAGCAAGTGCTGTCGCCGTTGAAGACTAC" \
    "GGATCACTGTAAGACCGGCAATAGATGAGCGCAACAGTCTACAATTAGACACCAGCATCA" \
    "CAAGAGGTAACAATCTCCACGCAGCCAATGCGGTCGTTTATTATCCGATGAGTGACGCCT" \
    "GGTCGACGAAAAACAAACCTAAACAATGTCACGTGCGGGTAAGTTCTCACTCGCCATTTA" \
    "CAGTAATGCCTATAACTAACAGACTTGCACCTGTTGAGGAACTATAGTGCACCCAACGTA" \
    "AGCTGATGAAGAGAATAATTCTGTTTTATTCATTACTGTTTTGCATCGCCGTCAGCCCAT" \
    "TTGGCTCCGGTACGACAGGGAGGAGC"

real_double_barcode_read = "ATCATTGTGCTAGTTCGGTTACATTATTGCTAAGGTTAAGACACTTTCTGCCTTTGCGGAGAAACAGCACCTGGTGCTGCCTGGCTCTGGAAGAAAGAATGGACTTTTAAACCTACTTGCCTGTCGCTCTATCTTCAGCGTCTGCCCGGGTGTTTAACCTGAGTTTGATCTGGCTCGGGACGAACTGGCAGCGTGCCTAATAATGCAAGTCGAGCGGACTTAAAAAGCTTACTTTTAAGTTGGCAACGGACAGGTGAAATGGCACGTAAGACGACCTGCCTGTAAGACTGGGATAACTTCCGGGGAAACCAGGAGCTAATACCGGATAATCCTTTCCTACTCATGTAGGAAAGCTGAGAAGACTTTACGCTGTCTACTTACAGATGGGCCCGCAGCGCATTAGCAGTTGGTGGGTAACGGCTCCTGGACGACGATCGCGGCATACGACCTGAAGGGTGATCAGCCACACTGGGACTGAGACGGCGAACTCCTGGGAGCAGCAGTAAGAATCTTCCGCAATGGACGAAGTCTGACAGACAACGCACCGCGTAATTGATGAAGGTTTTTTCGGATCGTAAACTCTGTTGTTAGGGAAGAACAAGTACGAGAGTAACTGCTCGTACCTTGAGCGGTACCTAACCAAGCGGCTAACTACGTGCCAGCAACCGCGGTAATACGTAAGGTAGCAGAAGCGTTGTCCGGGATATTGGGCGTAAAGCGCTTATGAGCGATCGCTAGTCTGATGTGAAAGCCCACAGCTCAACCGTGGAGGTCATTGAAGCAGGGGACTTAAGGTACAGAAGAAGAGTGGAATTCCACGATGTAGCGGTGAAATGCGTAGAGATGTGGGGAACACCAATTGGCGGGCAACTCTTTGGTCTGTAACTGACGCTGAGGCGGCGGCGGCGGCGCGAAAGCGTGGGGAGCAAACAGGTTAGATACCCTGGTAGTCCACGCCGTAAACGATGGTGCTAAGTGTTGCAGGGTTTCCGCCCTTTGGTGCTGCAACAAGCGCGACAAGCACTCGCCTGGGGGTACGGCCGCAAGGCTGAAACTCAAAGATTGACGGGGGCCGCACAAGCGGTGGGCATGTGTGGTTTTAATTCGGAAGCAACGCGAAGAACCTTACCAGGTCTTGACATCCTCGGCAATCCTAGAGATAGGTTCCCCGCTTCAGGAGGGACAGGATGACAGTGGTGCATGGTGTCGTCAGCTCGTGATAAGTGAGATGTTGGGTTAAGTCCCGCAACAGGCACAGCCCTTGATCTTAGTTGCCAGCATTTAGTTGGGCACTCTAAGTGACGCCGGTGACAAACCGGAGGAGGAGGTGGGGATGACGTCAAATCATCATGCCCCTATGACCTGGGCCACTTGCAGTACTGCAATGGATAATTTTCAAAGGGCAGCAAAACCGCGAGGTCGAGCCAATCCCACAAAATAACCATTCTCGAGTTCACGGATTGCGGCTGCAGCTGCCTACATGGGCTGGAATCCTTTAGTAATCGCAGCGGATGGCATGCCGCGGTGAATACGTTCCCGGGCCTTCATACGCACCGCCGTCACACCACGAGTTTATGGCTGCAGAAGTCGGTGGGGTAACCTTCACAGAACCAGCCACTAGGTGGGATAGATGATTGGGGTGGCAAGTCGTAACAAGGTAACGAGGTTAAACACCCAAGCAGACGCCGCAATATCAGCAACAGAAGGTTAAAAGTCCATTCTTCTTCCAGACAGGCAGCACCGATTGCTGTTCTCGCAAGGTAAGAAAGTAGTCTTAACCTTAGCCATGCATGC"


def test_scanner_find_best_adapter_template():
    qcat_config = config.get_default_config()

    best_adapter_template, best_adapter_end_position, best_adapter_score =\
        find_best_adapter_template(None, None, qcat_config)

    assert best_adapter_template == -1
    assert best_adapter_end_position == -1

    best_adapter_template, best_adapter_end_position, best_adapter_score = \
        find_best_adapter_template(get_adapter_by_name("RBK001"), read_bc3_exact, qcat_config)

    assert best_adapter_template == 0
    assert best_adapter_end_position == 101


def test_barcode_rapidkit():

    detector = BarcodeScannerEPI2ME(kit="RBK001")

    barcode_dict = detector.detect_barcode(read)
    assert barcode_dict['barcode'].name == "barcode02"
    barcode_dict = detector.detect_barcode(read_bc3_exact)
    assert barcode_dict['barcode'].name == "barcode03"
    barcode_dict = detector.detect_barcode(read_bc3)
    assert barcode_dict['barcode'].name == "barcode03"
    barcode_dict = detector.detect_barcode(real_bc03_porechop)
    assert barcode_dict['barcode'].name == "barcode03"
    barcode_dict = detector.detect_barcode(read_nobc)
    assert barcode_dict['barcode'] is None
    barcode_dict = detector.detect_barcode("")
    assert barcode_dict['barcode'] is None


# def test_scanner_detect_barcode_simple():
#
#     detector = scanner.BarcodeScannerSimple()
#
#     barcode_dict = detector.detect_barcode(read)
#     assert barcode_dict['barcode'].name == "barcode02"
#     barcode_dict = detector.detect_barcode(read_bc3_exact)
#     assert barcode_dict['barcode'].name == "barcode03"
#     barcode_dict = detector.detect_barcode(read_bc3)
#     assert barcode_dict['barcode'].name == "barcode03"
#     barcode_dict = detector.detect_barcode(real_bc03_porechop)
#     assert barcode_dict['barcode'].name == "barcode03"
#     barcode_dict = detector.detect_barcode(read_nobc)
#     assert barcode_dict['barcode'] is None
#     barcode_dict = detector.detect_barcode("")
#     assert barcode_dict['barcode'] is None
#
#     barcode_dict = detector.detect_barcode(read)
#     assert barcode_dict['barcode'].name == "barcode02"
#     barcode_dict = detector.detect_barcode(read_bc3_exact)
#     assert barcode_dict['barcode'].name == "barcode03"
#     barcode_dict = detector.detect_barcode(read_bc3)
#     assert barcode_dict['barcode'].name == "barcode03"
#     barcode_dict = detector.detect_barcode(real_bc03_porechop)
#     assert barcode_dict['barcode'].name == "barcode03"


def test_scanner_detect_barcode_generic():
    qcat_config = config.get_default_config()
    # config.min_quality = 0.1

    detector = BarcodeScannerEPI2ME()

    seq = extract_align_sequence(read, False, qcat_config.max_align_length)
    barcode_dict = detector.detect_barcode(seq, qcat_config=qcat_config)
    assert barcode_dict['barcode'].name == "barcode02"

    seq = extract_align_sequence(read_bc3_exact, False, qcat_config.max_align_length)
    barcode_dict = detector.detect_barcode(seq, qcat_config=qcat_config)
    assert barcode_dict['barcode'].name == "barcode03"

    seq = extract_align_sequence(read_bc3, False,
                                         qcat_config.max_align_length)
    barcode_dict = detector.detect_barcode(seq, qcat_config=qcat_config)
    assert barcode_dict['barcode'].name == "barcode03"

    seq = extract_align_sequence(real_bc03_porechop, False,
                                         qcat_config.max_align_length)
    barcode_dict = detector.detect_barcode(seq)
    assert barcode_dict['barcode'].name == "barcode03"

    seq = extract_align_sequence(read_nobc, False,
                                         qcat_config.max_align_length)

    barcode_dict = detector.detect_barcode(seq)
    assert barcode_dict['barcode'] is None

    seq = extract_align_sequence("", False,
                                         qcat_config.max_align_length)
    barcode_dict = detector.detect_barcode(seq)
    assert barcode_dict['barcode'] is None

    seq = extract_align_sequence(simulated_bc03, False,
                                         qcat_config.max_align_length)
    barcode_dict = detector.detect_barcode(seq)
    assert barcode_dict['barcode'].name == "barcode03"

    seq = extract_align_sequence(real_bc03_porechop, False,
                                         qcat_config.max_align_length)
    barcode_dict = detector.detect_barcode(seq, qcat_config=qcat_config)
    assert barcode_dict['barcode'].name == "barcode03"

    seq = extract_align_sequence(real_bc03_porechop, False,
                                         qcat_config.max_align_length)
    barcode_dict = detector.detect_barcode(seq, qcat_config=qcat_config)
    assert barcode_dict['barcode'].name == "barcode03"


# def test_config():
#
#     config_write = config.get_default_config()
#
#     config_write.match = 5
#     config_write.mismatch = -5
#     config_write.gap_open = -5
#     config_write.gap_extend = -5
#
#     config_write.write("test.ini")
#
#     config_read = config.get_default_config()
#     config_read.read("test.ini")
#
#     assert config_write.match == config_read.match
#     assert config_write.mismatch == config_read.mismatch
#     assert config_write.gap_open == config_read.gap_open
#     assert config_write.gap_extend == config_read.gap_extend
#     assert config_write.max_align_length == config_read.max_align_length
#     assert config_write.extracted_barcode_extension == config_read.extracted_barcode_extension
#     assert config_write.extracted_barcode_extension == config_read.extracted_barcode_extension


def test_utils():
    assert utils.qstring_to_phred(None) == []
    assert utils.qstring_to_phred("") == []
    assert utils.qstring_to_phred("IIII") == [40, 40, 40, 40]

    assert utils.mean_error_prob(None) == -1
    assert utils.mean_error_prob([]) == -1
    assert utils.mean_error_prob([40]) == 0.0001
    assert utils.mean_error_prob([40, 40, 40, 40]) == 0.0001
    assert utils.mean_error_prob([40, 20, 20, 40]) == 0.00505


def test_scanner_ectract_align_sequence():
    align_sequence = extract_align_sequence(None, True, 100)
    assert align_sequence == ""
    align_sequence = extract_align_sequence(None, False, 100)
    assert align_sequence == ""
    align_sequence = extract_align_sequence(None, True, 0)
    assert align_sequence == ""
    align_sequence = extract_align_sequence(None, False, 0)
    assert align_sequence == ""
    align_sequence = extract_align_sequence(None, True, -1)
    assert align_sequence == ""
    align_sequence = extract_align_sequence(None, False, -1)
    assert align_sequence == ""

    align_sequence = extract_align_sequence("", True, 100)
    assert align_sequence == ""
    align_sequence = extract_align_sequence("", False, 100)
    assert align_sequence == ""
    align_sequence = extract_align_sequence("", True, 0)
    assert align_sequence == ""
    align_sequence = extract_align_sequence("", False, 0)
    assert align_sequence == ""
    align_sequence = extract_align_sequence("", True, -1)
    assert align_sequence == ""
    align_sequence = extract_align_sequence("", False, -1)
    assert align_sequence == ""

    read_sequence = "AGTATTACTTCGTTCAGTTACGTATTGCTGTTTCATCTATCAGGAGGGAATGGAGTTTCGC"
    align_sequence = extract_align_sequence(read_sequence, False, 10)
    assert align_sequence== "AGTATTACTT"
    align_sequence = extract_align_sequence(read_sequence, True, 10)
    assert align_sequence == "GCGAAACTCC"
    align_sequence = extract_align_sequence(read_sequence, True, 0)
    assert align_sequence == read_sequence
    align_sequence = extract_align_sequence(read_sequence, False, 0)
    assert align_sequence == read_sequence
    align_sequence = extract_align_sequence(read_sequence, True, -1)
    assert align_sequence == read_sequence
    align_sequence = extract_align_sequence(read_sequence, False, -1)
    assert align_sequence == read_sequence


# def test_scanner_compute_adapter_error_prob():
#     assert scanner.compute_adapter_error_prob(None, 3, 3) == -1.0
#     assert scanner.compute_adapter_error_prob("", 3, 3) == -1.0
#     assert scanner.compute_adapter_error_prob(None, -3, -3) == -1.0
#     assert scanner.compute_adapter_error_prob("", -3, -3) == -1.0
#
#     assert scanner.compute_adapter_error_prob("III", 3, 3) == 0.0001
#     assert scanner.compute_adapter_error_prob("III", 300, 3) == -1.0
#     assert scanner.compute_adapter_error_prob("III", -300, 3) == -1.0
#     assert scanner.compute_adapter_error_prob("III", -300, -3) == -1.0


def test_adapter_layout():

    # test_layout = \
    #     barcodes.AdapterLayout("adapter_layout_pbc001_pbc096", 1,
    #                   "AAAAAAAAAT" + barcodes.barcode_spacer + "ATTTTTTTTT", 1,
    #                   adapters.get_barcodes_simple(),
    #
    #    "GGGGGGGGGC" + barcodes.barcode_spacer + "GCCCCCCCCC", 1,
    #                   barcodes.barcodes_standard_rev,
    #                   [barcodes.barcoding_kits["PBC001"]],
    #                   "Test layout")

    test_layout_5p = adapters.AdapterLayout(
                               "PBC001",
                               "AAAAAAAAAT" + barcode_spacer + "ATTTTTTTTT",
                               adapters.get_barcodes_simple(),
                               None,
                               "Test layout")
    test_layout_3p = adapters.AdapterLayout(
                               "PBC001",
                               "GGGGGGGGGC" + barcode_spacer + "GCCCCCCCCC",
                               adapters.get_barcodes_simple(),
                               None,
                               "Test layout")
    test_layout_double = adapters.AdapterLayout(
                                                "PBC001",
                                                "AAAAAAAAAT" +
                                                barcode_spacer +
                                                "ATTTTTTTTTGGGGGGGGGC" +
                                                barcode_spacer +
                                                "GCCCCCCCCC",
                                                adapters.get_barcodes_simple(),
                                                adapters.get_barcodes_simple(),
                                                "Test layout")

    assert test_layout_double.get_adapter_sequences() == ("AAAAAAAAAT" + barcode_spacer + "ATTTTTTTTTGGGGGGGGGC" + barcode_spacer + "GCCCCCCCCC")

    assert test_layout_double.get_adapter_length() == len("AAAAAAAAAT" + barcode_spacer + "ATTTTTTTTTGGGGGGGGGC" + barcode_spacer + "GCCCCCCCCC")

    assert test_layout_double.get_barcode_end(0) == (len("AAAAAAAAAT") + len(barcode_spacer) - 1)
    assert test_layout_double.get_barcode_end(1) == len("AAAAAAAAAT" + barcode_spacer + "ATTTTTTTTTGGGGGGGGGC" + barcode_spacer) - 1

    assert test_layout_double.get_barcode_length(0) == len(barcode_spacer)
    assert test_layout_double.get_barcode_length(1) == len(barcode_spacer)

    assert test_layout_double.get_barcode_set(0) == adapters.get_barcodes_simple()
    assert test_layout_double.get_barcode_set(1) == adapters.get_barcodes_simple()

    assert test_layout_double.get_downstream_context(2, 0) == "AT"
    assert test_layout_double.get_downstream_context(2, 1) == "GC"

    assert test_layout_double.get_upstream_context(2, 0) == "AT"
    assert test_layout_double.get_upstream_context(2, 1) == "GC"


def _parse_reads_info(comment):
    single_read_info = {}
    if comment:
        for pair in comment.split(" "):
            cols = pair.split('=')
            single_read_info[cols[0]] = cols[1]
    return single_read_info


def full_run(path):

    correct = 0
    incorrect = 0
    fn = 0
    fp = 0

    detector = scanner.factory()
    for names, comments, seqs, quals in cli.iter_fastx(path, True, 1):

        name, comment, sequence, quality = names[0], comments[0], seqs[0], quals[0]

        result = detector.detect_barcode(sequence)

        if result:
            trim5p = result['trim5p']
            trim3p = result['trim3p']

            trimmed_len = trim3p - trim5p

            # assert trimmed_len > 0
            assert len(sequence) >= trimmed_len
            assert trim5p < 200
            assert (len(sequence) - trim3p) < 200

        info = _parse_reads_info(comment)
        truebc = info['truebc']
        bc = "none"
        if result['barcode']:
            bc = str(result['barcode'].id)

        print(bc, truebc)
        if bc == truebc:
            correct += 1
        else:
            if bc == "none":
                fn += 1
            else:
                if truebc == "none":
                    fp += 1
                else:
                    incorrect += 1

    total = correct + incorrect + fn + fp

    incorrect_perc = (incorrect + fp) * 100.0 / total
    fn_perc = fn * 100.0 / total

    return incorrect_perc, fn_perc


def test_full_run_rad003():

    path = "qcat/test/data/nobarcode_1k.fastq"
    incorrect_perc, fn_perc = full_run(path)

    print("Incorrect: {}, False negative: {}".format(incorrect_perc, fn_perc))
    assert incorrect_perc <= 0.08
    assert fn_perc == 0.0


def test_full_run_nbd103():

    path = "qcat/test/data/nbd103.fastq"
    incorrect_perc, fn_perc = full_run(path)

    print("Incorrect: {}, False negative: {}".format(incorrect_perc, fn_perc))
    assert incorrect_perc <= 0.00
    assert fn_perc == 0.0


def test_full_run_pbk004():

    path = "qcat/test/data/pbk004.fastq"
    incorrect_perc, fn_perc = full_run(path)

    print("Incorrect: {}, False negative: {}".format(incorrect_perc, fn_perc))
    assert incorrect_perc <= 0.00
    assert fn_perc == 0.0


def test_full_run_rab204():

    path = "qcat/test/data/rab204.fastq"
    incorrect_perc, fn_perc = full_run(path)

    print("Incorrect: {}, False negative: {}".format(incorrect_perc, fn_perc))
    assert incorrect_perc <= 0.00
    assert fn_perc == 0.0


def test_full_run_rbk004():

    path = "qcat/test/data/rbk004.fastq"
    incorrect_perc, fn_perc = full_run(path)

    print("Incorrect: {}, False negative: {}".format(incorrect_perc, fn_perc))
    assert incorrect_perc <= 0.00
    assert fn_perc == 0.0


# def test_full_run_nbd114():
#
#     path = "qcat/test/data/nbd_104_114_4k.fastq"
#     incorrect_perc, fn_perc = full_run(path)
#
#     print("Incorrect: {}, False negative: {}".format(incorrect_perc, fn_perc))
#     assert incorrect_perc <= 0.7  # 1.1
#     assert fn_perc <= 11.0


def test_trimming():

    fp = 0
    tn = 0

    path = "qcat/test/data/barcode_1k.fastq"
    i = 0

    demultiplexer = BarcodeScannerEPI2ME()

    for names, comments, seqs, quals in cli.iter_fastx(path, True, 1):

        name, comment, sequence, quality = (
        names[0], comments[0], seqs[0], quals[0])

        result = demultiplexer.detect_barcode(sequence)

        assert result["trim3p"] > result["trim5p"]
        if result and result['barcode']:
            trim_5p = result["trim5p"]
            trim_3p = result["trim3p"]
            sequence = sequence[trim_5p:trim_3p]
            assert len(sequence) == result["trim3p"] - result["trim5p"]
            result = demultiplexer.detect_barcode(sequence)
            if result['barcode']:
                fp += 1
            else:
                tn += 1
        else:
            pass
        i += 1

    fp = fp / i * 100.0

    assert fp <= 0.2 # 2.0

