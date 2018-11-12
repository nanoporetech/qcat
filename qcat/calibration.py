# import pkg_resources
# # import mappy
#
# RESOURCE_PACKAGE = __name__
# LAMBDA_DNA_MM2_IDX = '/'.join(('resources', 'lambda_3.6kb.fasta.minimap2.idx'))
# LAMBDA_DNA_MM2_IDX = pkg_resources.resource_filename(__name__, LAMBDA_DNA_MM2_IDX)
#
#
# def detect_calibration_strand(read_sequence, is_dna):
#     """
#     Checks whether the read_sequence stems from the lambda calibration strand
#     :param read_sequence: sequence of the read
#     :param is_dna: DNA or RNA read
#     :return: dict
#     """
#
#     if is_dna:
#         if detect_calibration_strand.minimap2_lambda_dna is None:
#             detect_calibration_strand.minimap2_lambda_dna = \
#                 mappy.Aligner(LAMBDA_DNA_MM2_IDX)
#         minimap2 = detect_calibration_strand.minimap2_lambda_dna
#     else:
#         raise NotImplementedError("RNA calibration strand detection not"
#                                   "implemented yet.")
#
#     mapped = False
#     r_len = len(read_sequence)
#     q_st = -1
#     q_en = -1
#     hits = minimap2.map(read_sequence)
#     for hit in hits:
#         if hit.is_primary:
#             q_st = hit.q_st
#             q_en = hit.q_en
#             qry_cov = (hit.q_en - hit.q_st) * 100.0 / r_len
#             if qry_cov > 90.0:
#                 mapped = True
#
#     return_dict = {"calibration_strand": mapped,
#                    "qry_start": q_st,
#                    "qry_end": q_en,
#                    "exit_status": 0
#                   }
#     return return_dict
#
#
# detect_calibration_strand.minimap2_lambda_dna = None
