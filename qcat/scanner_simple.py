import logging
import os

from qcat import config
from qcat.adapters import get_barcodes_simple, get_barcodes_from_fastq
from qcat.scanner_base import BarcodeScanner, find_highest_scoring_barcode, \
    build_return_dict, empty_return_dict


class BarcodeScannerSimple(BarcodeScanner):

    def __init__(self, min_quality=None, kit_folder=None, kit=None, enable_filter_barcodes=False, scan_middle_adapter=False, threads=1):

        if min_quality is None:
            min_quality = 60

        if threads != 1:
            logging.warning("Multi threading is not yet supported in "
                            "simple mode. Falling back to using a "
                            "single thread.")

        super(BarcodeScannerSimple, self).__init__(min_quality,
                                                   None,
                                                   kit_folder=kit_folder,
                                                   enable_filter_barcodes=enable_filter_barcodes,
                                                   scan_middle_adapter=scan_middle_adapter
                                                   )
        if os.path.isfile(kit) and os.path.exists(kit):
            self.barcodes = get_barcodes_from_fastq(kit)
        else:
            self.barcodes = get_barcodes_simple(kit)

    @staticmethod
    def get_name():
        return "simple"

    def barcode_count(self):
        if self.barcodes:
            return len(self.barcodes) + 1
        else:
            return super(BarcodeScannerSimple, self).barcode_count()

    def scan(self,
             read_sequence,
             read_qualities,
             bc_adapter_templates,
             nobc_adapter_templates,
             qcat_config=config.qcatConfig()):
        """
            Aligns all passed barcodes to the read sequence and chooses the one
            with the highest alignment score. Simplest version of barcode detection
            that is independent of the sequencing kit but has higher FN, FP and
            incorrect rate for barcode assignments. Currently does not implement any
            identity thresholds. Therefor, filtering by quality score is required.

            :param read_sequence: Sequence of the read
            :type read_sequence: str
            :param read_qualities: Base qualities of the read in Sanger encoding
            :type read_qualities: str
            :param barcode_set: List of possible Barcode tuples
            :type barcode_set: List
            :param qcat_config: qcatConfig object
            :type qcat_config: qcatConfig
            :return: see build_return_dict
            :rtype: Dictionary
            """

        exit_status = 0

        best_barcode, best_barcode_q_score, identity, barcode_end = \
            find_highest_scoring_barcode(barcode_region_read=read_sequence,
                                         barcode_set=self.barcodes,
                                         qcat_config=qcat_config,
                                         compute_identity=True)

        # barcode_err_probe = 0.0
        # if best_barcode:
        #     barcode_err_probe = compute_adapter_error_prob(read_qualities,
        #                                                    barcode_end,
        #                                                    len(
        #                                                        best_barcode.sequence))

        # If identity is too low or error prob too high, report as unclassified
        if identity < self.min_quality:
            return empty_return_dict()

        return build_return_dict(
            best_barcode=best_barcode,
            best_barcode_score=best_barcode_q_score,
            best_adapter=None,
            best_adapter_end=barcode_end,
            exit_status=exit_status
        )