import logging

from qcat import config
from qcat.scanner import BarcodeScanner
from qcat.scanner_base import find_best_adapter_template, \
    extract_barcode_region, find_highest_scoring_barcode, build_return_dict


class BarcodeScannerEPI2ME(BarcodeScanner):

    def __init__(self, min_quality=None, kit_folder=None, kit=None, enable_filter_barcodes=False, scan_middle_adapter=False, threads=1):

        if min_quality is None:
            min_quality = 58

        if threads != 1:
            logging.warning("Multi threading is not yet supported in "
                            "epi2mme mode. Falling back to using a "
                            "single thread.")

        super(BarcodeScannerEPI2ME, self).__init__(min_quality,
                                                   kit,
                                                   kit_folder=kit_folder,
                                                   enable_filter_barcodes=enable_filter_barcodes,
                                                   scan_middle_adapter=scan_middle_adapter,
                                                   )
        self.barcodes = None

    @staticmethod
    def get_name():
        return "epi2me"

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
        if best_adapter_score > 90.0 or best_adapter_template.is_doulbe_barcode():
            barcode_region_read = extract_barcode_region(
                read_sequence=read_sequence,
                adapter_template=best_adapter_template,
                barcode_set_index=barcode_set_index,
                alignment_stop_ref=aligned_adapter_end,
                qcat_config=qcat_config)
        else:
            barcode_region_read = read_sequence[:qcat_config.max_align_length]

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

        # If double barcode adapter
        barcode_set_index = 1
        if best_adapter_template.get_barcode_set(barcode_set_index):

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
        else:
            logging.debug("Adapter type does not have second barcode")

        aligned_adapter_end = min(aligned_adapter_end +
                                  best_adapter_template.trim_offset,
                                  len(read_sequence))
        return build_return_dict(
            best_barcode=best_barcode,
            best_barcode_score=best_barcode_q_score,
            best_adapter=best_adapter_template,
            best_adapter_end=aligned_adapter_end,
            exit_status=exit_status
        )
