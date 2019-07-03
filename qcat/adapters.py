from builtins import AttributeError

import glob
import os
import logging

import pkg_resources
import yaml
from collections import namedtuple

from Bio.SeqIO.FastaIO import SimpleFastaParser

from qcat.layout import AdapterLayout

Barcode = namedtuple("Barcode", "name id sequence fwd_strand")


def barcode2yaml(bc):
    barcode = {}
    barcode['name'] = bc.name
    barcode['id'] = bc.id
    barcode['sequence'] = bc.sequence
    barcode['fwd_strand'] = bc.fwd_strand
    return barcode


def barcodes2yaml(barcodes):
    barcode_set = []
    for bc in barcodes or []:
        barcode_set.append(barcode2yaml(bc))
    return barcode_set


def yaml2adapter(filename):
    with open(filename, 'r') as stream:
        try:
            data = yaml.load(stream)
            return data
        except yaml.YAMLError as exc:
            print(exc)


def adapter2yaml(adapter):
    return {"kit": adapter.kit,
            "auto_detect": adapter.auto_detect,
            "model": adapter.model,
            "description": adapter.description,
            "sequence": adapter.get_adapter_sequences(),
            "barcode_set_1": barcodes2yaml(adapter.get_barcode_set(0)),
            "barcode_set_2": barcodes2yaml(adapter.get_barcode_set(1))}


def read_barcode(data):
    if not data:
        return None
    return Barcode(data['name'], data['id'], data.get('sequence', None), data.get('fwd_strand', None))


def read_barcode_set(data):
    if not data:
        return None

    barcodes = []
    for barcode in data:
        barcodes.append(read_barcode(barcode))

    return barcodes


RESOURCE_PACKAGE = __name__
KIT_FOLDER = '/'.join(('resources', 'kits'))
KIT_FOLDER = pkg_resources.resource_filename(__name__, KIT_FOLDER)


def read_adapter_layout(filename):
    with open(filename, 'r') as stream:
        try:
            try:
                # Python 3
                data = yaml.load(stream, Loader=yaml.FullLoader)
            except AttributeError:
                # Python 2
                data = yaml.load(stream)

            if not data.get('active', True):
                return None

            model = None
            model_len = None
            if "model" in data and data['model']:
                model = data['model'].get('file', None)
                model_len = data['model'].get('length', None)

            return AdapterLayout(
                kit=data.get('kit', ""),
                sequence=data.get('sequence', ""),
                barcode_set_1=read_barcode_set(data.get('barcode_set_1', None)),
                barcode_set_2=read_barcode_set(data.get('barcode_set_2', None)),
                description=data.get('description', ""),
                auto_detect=data.get('auto_detect', False),
                trim_offset=data.get("trim_offset", 0),
                model=model,
                model_len=model_len)
        except yaml.YAMLError as exc:
            print(exc)


def get_barcodes_from_fastq(reads_fa):
    barcode_id = 1
    barcodes = []
    with open(reads_fa) as f:
        for title, seq in SimpleFastaParser(f):
            barcodes.append(read_barcode({'name': title, 'id': barcode_id, 'sequence': seq}))
            barcode_id += 1

    if len(barcodes) <= 0:
        logging.error("Couldn't find barcodes in {}".format(reads_fa))
    return barcodes


def get_barcodes_simple(kit="standard", filename=None):
    if not filename or not os.path.isfile(filename):
        filename = os.path.join(KIT_FOLDER, "simple_{}.yml".format(kit))

    with open(filename, 'r') as stream:
        try:
            try:
                # Python 3
                data = yaml.load(stream, Loader=yaml.FullLoader)
            except AttributeError:
                # Python 2
                data = yaml.load(stream)
            return read_barcode_set(data.get('barcode_set_1', []))
        except yaml.YAMLError as exc:
            print(exc)


def populate_adapter_layouts(folder=None):
    if not folder:
        folder = KIT_FOLDER

    if os.path.exists(folder):
        if os.path.isdir(folder):
            filenames = glob.glob(os.path.join(folder, "*.yml"))
        else:
            filenames = [folder]
    else:
        logging.warn("{} not found. Using default adapter sequences.".format(folder))
        filenames = [KIT_FOLDER]

    adapters = []
    for filename in filenames:
        # print(filename)
        layout = read_adapter_layout(filename)
        # print(layout)
        if not layout:
            continue
        if layout.auto_detect:
            adapters.append(layout)
        else:
            adapters.append(layout)
    return adapters
