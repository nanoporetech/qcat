import logging

from qcat import config
from qcat.adapters import Barcode
from qcat.scanner import BarcodeScanner
from qcat.scanner_base import find_best_adapter_template, \
    extract_barcode_region, find_highest_scoring_barcode, build_return_dict, \
    empty_return_dict


class BarcodeScannerDual(BarcodeScanner):

    def __init__(self, min_quality=None, kit_folder=None, kit=None, enable_filter_barcodes=False, scan_middle_adapter=False, threads=1):

        if min_quality is None:
            min_quality = 60

        if threads != 1:
            logging.warning("Multi threading is not yet supported in "
                            "dual mode. Falling back to using a "
                            "single thread.")

        super(BarcodeScannerDual, self).__init__(min_quality,
                                                 "dual",
                                                 kit_folder=kit_folder,
                                                 enable_filter_barcodes=enable_filter_barcodes,
                                                 scan_middle_adapter=scan_middle_adapter
                                                 )
        self.barcodes = None

    @staticmethod
    def get_name():
        return "dual"

    def scan(self,
             read_sequence,
             read_qualities,
             bc_adapter_templates,
             nobc_adapter_templates,
             qcat_config=config.qcatConfig()):
        """
            Detects sequencing adapter and identifies best matching barcode

            :param read_sequence: Read sequence containing adapter
            :type read_sequence: str
            :param read_qualities: Base qualities of the read in Sanger encoding
            :type read_qualities: str
            :param adapter_templates: List of AdapterLayout objects
            with Ns
            :type adapter_templates: List
            :param qcat_config: qcatConfig object
            :type qcat_config: qcatConfig
            :return: see build_return_dict
            :rtype: Dictionary
            """
        exit_status = 0

        # Finding best barcoded adapters
        ret = find_best_adapter_template(adapter_templates=bc_adapter_templates,
                                         read_sequence=read_sequence,
                                         qcat_config=qcat_config)

        best_adapter_template_index, aligned_adapter_end, best_adapter_score = ret

        # if barcoded adapter found
        best_adapter_template = bc_adapter_templates[
            best_adapter_template_index]

        # Detect barcode
        barcode_set_index = 0

        # For high quality adapter alignments just use the barcode region,
        # for low quality compare full adapter to the barcodes
        # Do not scann full adapter is it contains double barcoding
        # best_adapter_score > 90.0 or

        barcode_region_read = extract_barcode_region(
            read_sequence=read_sequence,
            adapter_template=best_adapter_template,
            barcode_set_index=barcode_set_index,
            alignment_stop_ref=aligned_adapter_end,
            qcat_config=qcat_config)

        # First barcode
        up_context = best_adapter_template.get_upstream_context(
            qcat_config.barcode_context_length, barcode_set_index)
        down_context = best_adapter_template.get_downstream_context(
            qcat_config.barcode_context_length, barcode_set_index)
        barcode_set = best_adapter_template.get_barcode_set(
            barcode_set_index)
        if self.barcodes:
            barcode_set = self.barcodes

        # Find best barcode
        best_barcode, best_barcode_q_score, best_barcode_score, barcode_end = \
            find_highest_scoring_barcode(
                barcode_region_read=barcode_region_read,
                barcode_set=barcode_set,
                upstream_context=up_context,
                downstream_context=down_context,
                qcat_config=qcat_config)

        # Detect second barcode
        barcode_set_index = 1
        barcode_region_read = extract_barcode_region(
            read_sequence=read_sequence,
            adapter_template=best_adapter_template,
            barcode_set_index=barcode_set_index,
            alignment_stop_ref=aligned_adapter_end,
            qcat_config=qcat_config)

        # First barcode
        up_context = best_adapter_template.get_upstream_context(
            qcat_config.barcode_context_length, barcode_set_index)
        down_context = best_adapter_template.get_downstream_context(
            qcat_config.barcode_context_length, barcode_set_index)
        barcode_set = best_adapter_template.get_barcode_set(
            barcode_set_index)
        if self.barcodes:
            barcode_set = self.barcodes

        # Find best barcode
        best_barcode_2, best_barcode_q_score_2, best_barcode_score_2, barcode_end_2 = \
            find_highest_scoring_barcode(
                barcode_region_read=barcode_region_read,
                barcode_set=barcode_set,
                upstream_context=up_context,
                downstream_context=down_context,
                qcat_config=qcat_config)

        if best_barcode and best_barcode_2:
            dual_barcode = Barcode("barcode{:02d}/{:02d}".format(best_barcode.id, best_barcode_2.id),
                                   "{}/{}".format(best_barcode.id, best_barcode_2.id),
                                   None,
                                   True
                                   )

            return build_return_dict(
                best_barcode=dual_barcode,
                best_barcode_score=min(best_barcode_q_score, best_barcode_q_score_2),
                best_adapter=best_adapter_template,
                best_adapter_end=aligned_adapter_end,
                exit_status=exit_status
            )
        else:
            return empty_return_dict()
