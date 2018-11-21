from __future__ import print_function

import sys
import logging

from qcat import config
from qcat.adapters import read_barcode
from qcat.layout import AdapterLayout
from qcat.scanner import BarcodeScanner

try:
    from pyguppy import barcoding
    guppy_import_failed = False
except ImportError as e:
    # logging.warning("Could not load pyguppy, skipping.")
    guppy_import_failed = True
    pass

from qcat.scanner_base import build_return_dict


class BarcodeScannerGuppy(BarcodeScanner):

    # This is not very nice. Should be part of the config files.
    KIT_MAPPING = {'NBD104/NBD114': 'NB',
                   'PBC096': 'BC',
                   'PBC001': 'BC',
                   'RBK001': 'RBK',
                   'RBK004': 'RBK4',
                   'PBK004/LWB001': 'LWB',
                   'RPB004/RLB001': 'RLB',
                   'RAB204': 'RAB'}

    KIT_MAPPING_REV = {'NB':   'NBD104/NBD114',
                       'BC':   'PBC096',
                       'RBK':  'RBK001',
                       'RBK4': 'RBK004',
                       'LBW':  'PBK004/LWB001',
                       'RLB':  'RPB004/RLB001',
                       'RAB':  'RAB204'}

    def __init__(self, min_quality=None, kit_folder=None, kit=None,
                 enable_filter_barcodes=False, scan_middle_adapter=False,
                 threads=1):

        if min_quality is None:
            min_quality = 60

        if scan_middle_adapter:
            logging.warning("Guppy/Albacore do not support scanning for "
                            "middle adapters. Ignoring parameter.")

        if enable_filter_barcodes:
            logging.warning("Guppy/Albacore do not support filtering "
                            "barcodes in batch mode. Ignoring parameter.")

        super(BarcodeScannerGuppy, self).__init__(min_quality,
                                                  kit,
                                                  kit_folder=kit_folder)

        # Determine if user specified a kit
        guppy_kits = set()
        for layout in self.layouts:
            guppy_kits.update([layout.kit])

        kit_name = None
        if len(guppy_kits) == 1:
            kit_name = BarcodeScannerGuppy.KIT_MAPPING.get(guppy_kits.pop())

        if kit_name:
            # print("Kit: {}".format(kit_name), file=sys.stderr)
            self.barcoder = barcoding.Barcoder(kit_name=kit_name,
                                               num_threads=threads,
                                               min_quality=min_quality)
        else:
            # print("Auto detecting kit", file=sys.stderr)
            self.barcoder = barcoding.Barcoder(num_threads=threads,
                                               min_quality=min_quality)

    @staticmethod
    def get_name():
        return "guppy"

    @staticmethod
    def convert_guppy_result(guppy_result):

        guppy_barcode = guppy_result['barcode']
        if guppy_barcode['name'] == 'unclassified':
            guppy_barcode = None
            guppy_adapter = None
        else:
            guppy_barcode = read_barcode(guppy_barcode)
            guppy_adapter = AdapterLayout(
                kit=BarcodeScannerGuppy.KIT_MAPPING_REV.get(guppy_result['kit'],
                                                            guppy_result['kit']),
                sequence="ATGC",
                barcode_set_1=None,
                barcode_set_2=None,
                description="")

        result = build_return_dict(
            best_barcode=guppy_barcode,
            best_barcode_score=guppy_result['score'],
            best_adapter=guppy_adapter,
            best_adapter_end=90,
            trim5p=guppy_result['trim5p'],
            trim3p=guppy_result['trim3p'],
            exit_status=guppy_result['exit_status'],
        )

        return result

    def detect_barcode(self,
                       read_sequence,
                       read_qualities=None,
                       qcat_config=config.qcatConfig()):

        guppy_result = self.barcoder.detect_barcode(("", read_sequence))[0]

        return self.convert_guppy_result(guppy_result)

    def detect_barcode_batch(self, read_sequences, read_qualities=[None],
                             qcat_config=config.qcatConfig()):

        guppy_results = self.barcoder.detect_barcode_batch([("", seq, "") for seq in read_sequences])

        results = []
        for guppy_result in guppy_results:
            results.append(self.convert_guppy_result(guppy_result))

        return results
