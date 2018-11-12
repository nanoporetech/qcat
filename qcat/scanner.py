from __future__ import print_function

import logging

from qcat import adapters
from qcat.scanner_base import BarcodeScanner
try:
    from qcat.scanner_guppy import BarcodeScannerGuppy, guppy_import_failed
    use_pyguppy = not guppy_import_failed
except ImportError as e:
    use_pyguppy = False

from qcat.scanner_epi2me import BarcodeScannerEPI2ME
from qcat.scanner_simple import BarcodeScannerSimple
from qcat.scanner_dual import BarcodeScannerDual


def get_adapter_by_name(kit, kit_folder=None):
    """
    Get Adapter object for Kit by name

    :param kit: Name of the Kit as given in the config file
    :param kit_folder: Folder containing config files for Adapters/Kits
    :return: List of Adapter objects with given name
    """
    adapter_list = []
    for adapter in adapters.populate_adapter_layouts(kit_folder):
        if adapter.kit == kit:
            adapter_list.append(adapter)
    return adapter_list


def get_modes():
    """
    Get all available algorithms for demultiplexing.
    Currently: qcat, brill, guppy, (simple)

    :return: List of available qcat modes
    """
    modes = []
    for subclass in BarcodeScanner.__subclasses__():
        name = subclass.get_name()
        if name != "brill" or use_brill or (name == "guppy" and use_pyguppy):
            modes.append(name)
    return modes


def get_kits_info(kit_folder=None):
    """
    Get name and description of all supported kits

    :param kit_folder: Folder containing config files for Adapters/Kits
    :return: Folder containing config files for Adapters/Kits
    """
    names = {'Auto': "Auto detect kit"}
    for layout in adapters.populate_adapter_layouts(kit_folder):
        if layout.kit not in names.keys():
            names[layout.kit] = layout.description

    return names


def get_kits(kit_folder=None):
    """
    Get names (as given in the config files) for all available kits.

    :param kit_folder: Folder containing config files for Adapters/Kits
    :return: List of kit names
    """
    names = ['Auto']
    for layout in adapters.populate_adapter_layouts(kit_folder):
        if layout.kit not in names:
            names.append(layout.kit)

    return names


def factory(mode="epi2me", min_quality=None, kit=None, kit_folder=None,
            enable_filter_barcodes=False, scan_middle_adapter=False, threads=1):
    """
    Create a BarcodeScanner object

    :param mode: Algorithm to use (qcat, brill, guppy, ...)
    :param min_quality: Minimum quality to call a barcode
    :param kit: Name of the kit
    :param kit_folder: Folder containing config files for Adapters/Kits
    :param enable_filter_barcodes: Remove barcodes that occur in small numbers
    when running in batch mode
    :param scan_middle_adapter: Scan full read for adapters
    :return: BarcodeScanner object
    """

    if mode == "guppy" and not use_pyguppy:
        logging.warning(
            "Demultiplexing mode {} currently not supported in your "
            "environment. Falling back to epi2me.".format(
                mode))
        mode = "epi2me"


    for subclass in BarcodeScanner.__subclasses__():
        if mode == subclass.get_name():
            return subclass(min_quality=min_quality,
                            kit_folder=kit_folder,
                            kit=kit,
                            enable_filter_barcodes=enable_filter_barcodes,
                            scan_middle_adapter=scan_middle_adapter,
                            threads=threads
                            )

    raise RuntimeError("Invalid demultiplexing mode: {}".format(mode))


# def detect_calibration_strand(read_sequence, is_dna):
#     """
#     Legacy function for running calibration strand detection
#
#     :param read_sequence: Sequence of the read
#     :param is_dna: DNA or RNA
#     :return: Calibration strand detected true/false
#     """
#     return calibration.detect_calibration_strand(read_sequence, is_dna)