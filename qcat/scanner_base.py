from __future__ import print_function

import logging
import sys

try:
    import parasail
except ImportError as e:
    logging.error("Could not load parasail library. Please try reinstalling "
                  "parasail using pip.")
    sys.exit(1)

import operator

from qcat import adapters, calibration
from qcat import config
from qcat.utils import revcomp


if parasail.can_use_sse2():
    parasail_sg_stat = parasail.sg_stats_striped_32
    parasail_sg = parasail.sg_striped_32
else:
    logging.warning("Warning: SSE not supported, falling back to standard alignment")
    parasail_sg_stat = parasail.sg_stats
    parasail_sg = parasail.sg


def extract_barcode_region(read_sequence, adapter_template, barcode_set_index,
                           alignment_stop_ref, qcat_config):
    """
    Extracts the region from the read that is expected to contain
    the barcode sequence.
    :param read_sequence: Sequence of the read
    :type read_sequence: str
    :param adapter_template: Adapter layout with masked barcode region
    :type adapter_template: AdapterLayout
    :param barcode_set_index:
    :param alignment_stop_ref: Position on the read that the last bp of the
    adapter template is aligned to
    :type alignment_stop_ref: int
    :param qcat_config: qcatConfig object
    :type qcat_config: qcatConfig
    :return: Part of the read that will be used to identify the barcode
    :rtype str
    """
    adapter_length = adapter_template.get_adapter_length()
    barcode_end = adapter_template.get_barcode_end(barcode_set_index)
    barcode_length = adapter_template.get_barcode_length(barcode_set_index)

    barcode_end_ref = alignment_stop_ref - (adapter_length - barcode_end) + 1
    barcode_start_ref = barcode_end_ref - barcode_length

    # Extend region to allow for indels
    barcode_start_ref -= min(qcat_config.extracted_barcode_extension,
                             barcode_start_ref)
    barcode_end_ref += min(qcat_config.extracted_barcode_extension,
                           len(read_sequence) - barcode_end_ref)

    return read_sequence[barcode_start_ref:barcode_end_ref + 1]


def find_highest_scoring_barcode(barcode_region_read,
                                 barcode_set,
                                 qcat_config,
                                 upstream_context="",
                                 downstream_context="",
                                 compute_identity=False):
    """
    Aligns all the barcodes from barcode_set to the barcode region that
    was extracted from the read and chooses the barcode with the best alignment
    score. Also computes a quality score for the barcode.
    Quality score = (best - secondbest) / best
    Where best is the highest alignment score found in all barcode
    alignments and secondbest is the second best score. Range is 0-1,
    a reasonable cutoff is between 0.1 and 0.2.

    :param barcode_region_read: Sequence extracted from the read that is
    supposed to contain the barcode
    :type barcode_region_read: str
    :param barcode_set: List of possible Barcode tuples
    :type barcode_set: List
    :param upstream_context: Sequence upstream of barcode in adapter_template
    :type upstream_context: str
    :param downstream_context: Sequence downstream of barcode in
    adapter_template
    :type downstream_context: str
    :param qcat_config: qcatConfig object
    :type qcat_config: qcatConfig
    :return: Best barcode and quality score
    :rtype Barcode, float
    """
    max_score = None
    max_identity = 0.0
    max_end = -1
    second_best_score = None
    second_best_identity = 0.0
    max_barcode = None
    q_score = 0

    if not barcode_region_read:
        return max_barcode, q_score, max_identity, max_end

    for barcode in barcode_set:

        if compute_identity:
            align = parasail_sg_stat
        else:
            align = parasail_sg

        aligned_barcode = align(s1=barcode_region_read,
                                s2=upstream_context +
                                   barcode.sequence +
                                   downstream_context,
                                open=1,
                                extend=1,
                                matrix=qcat_config.matrix_barcode)

        score = aligned_barcode.score * 100.0 / (1.0 * len(upstream_context + barcode.sequence + downstream_context))

        identity = 0.0
        if compute_identity:
            identity = float(aligned_barcode.matches) / len(barcode.sequence)

        if not max_score or max_score < score:
            second_best_score = max_score
            second_best_identity = max_identity
            max_score = score
            max_barcode = barcode
            max_identity = identity
            max_end = aligned_barcode.end_query
        elif not second_best_score or second_best_score < score:
            second_best_score = score
            second_best_identity = identity

    # print(max_score, file=sys.stderr)
    # if max_score > 0:
    #     q_score = int(float(max_score - second_best_score) / float(max_score) * 100.0)
    q_score = max_score

    return max_barcode, q_score, max_score, max_end


def align_adapter_identity(adapter_sequence, adapter_length, read_sequence,
                           barcode_length, qcat_config):
    """
    Aligns a single adapter template to the read an computes the
    identity for the alignment

    :param adapter_sequence: Sequence of adapter template
    :type adapter_sequence: str
    :param adapter_length: Length of adapter sequence
    :type adapter_length: int
    :param read_sequence: Read sequence
    :type read_sequence: str
    :param barcode_length: Length of the barcode. Number of
    bp will be ignored for identity computation
    :type barcode_length: int
    :param qcat_config: qcatConfig object
    :type qcat_config: qcatConfig
    :return: Parasail alignment object (None if partial match),
    alignment identity (0.0 if partial match)
    :rtype Result, int
    """
    if not read_sequence or not adapter_sequence:
        return None, 0.0

    aligned_adapter = parasail_sg_stat(s1=read_sequence,
                                       s2=adapter_sequence,
                                       open=qcat_config.gap_open,
                                       extend=qcat_config.gap_extend,
                                       matrix=qcat_config.matrix)

    # Check whether the whole adapter is aligned to the read or only a suffix
    partial_match = aligned_adapter.length < (adapter_length * 0.85)
    # partial_match = False
    if partial_match:
        # Currently partial matches between the read and the reference
        # are discarded
        aligned_adapter = None
        adapter_identity = 0.0
    else:
        # Compute alignment identity but do not take into account
        # the barcode region which is masked by Ns
        adapter_identity = float(aligned_adapter.matches) / \
                           float(aligned_adapter.length - barcode_length)

    return aligned_adapter, adapter_identity


def align_adapter(adapter_sequence, read_sequence, qcat_config):
    """
    Aligns a single adapter template to the read an computes the
    identity for the alignment

    :param adapter_sequence: Sequence of adapter template
    :type adapter_sequence: str
    :param adapter_length: Length of adapter sequence
    :type adapter_length: int
    :param read_sequence: Read sequence
    :type read_sequence: str
    :param barcode_length: Length of the barcode. Number of
    bp will be ignored for identity computation
    :type barcode_length: int
    :param qcat_config: qcatConfig object
    :type qcat_config: qcatConfig
    :return: Parasail alignment object (None if partial match),
    alignment identity (0.0 if partial match)
    :rtype Result, int
    """
    if not read_sequence or not adapter_sequence:
        return None, 0.0

    aligned_adapter = parasail_sg(s1=read_sequence,
                                  s2=adapter_sequence,
                                  open=qcat_config.gap_open,
                                  extend=qcat_config.gap_extend,
                                  matrix=qcat_config.matrix)

    return aligned_adapter, 0.0


def extract_align_sequence(read_sequence, rev_comp, length):
    """
    Extracts the part of the read that will be scanned for an adapter

    :param read_sequence: str
    :param length: Length of sequence that should be scanned. Starting at the 5'
    end of the read. If 0, full read is returned
    :type length: int
    :return: part read sequence
    :rtype str
    """
    if read_sequence:
        align_sequence = read_sequence
    else:
        align_sequence = ""
    if length > 0:
        if not rev_comp:
            align_sequence = align_sequence[:length]
        else:
            align_sequence = revcomp(align_sequence[-length:])

    return align_sequence


def compute_adapter_identity(adapter_template, read_sequence, qcat_config):
    aligned_adapter, adapter_identity = \
        align_adapter_identity(adapter_template.get_adapter_sequences(),
                               adapter_template.get_adapter_length(),
                               read_sequence,
                               adapter_template.get_barcode_length(
                                   0) + adapter_template.get_barcode_length(1),
                               qcat_config)
    return adapter_identity


def eval_adapter_template(adapter_template, read_sequence,
                          qcat_config, identity=True):
    """
    Extracts the adapter sequence with masked barcode from the adapter_template
    and aligns (semi-global) it to the 5' end, 3' end or full read.

    :param adapter_template: Adapter layout object
    :type adapter_template: AdapterLayout
    :param read_sequence: Read sequence
    :type read_sequence: str
    :param qcat_config: qcatConfig object
    :type qcat_config: qcatConfig
    :return: Position the last bp of the adapter is aligned to in the read,
    Identity of adapter alignment, Score of adapter alignment
    :rtype: int, float, float
    """
    adapter_end_position = -1
    adapter_score = -1

    # Align adapter to extracted sequence
    if identity:
        ret = align_adapter_identity(adapter_template.get_adapter_sequences(),
                                     adapter_template.get_adapter_length(),
                                     read_sequence,
                                     adapter_template.get_barcode_length(0) +
                                     adapter_template.get_barcode_length(1),
                                     qcat_config)
        aligned_adapter, adapter_identity = ret
    else:
        ret = align_adapter(adapter_template.get_adapter_sequences(),
                            read_sequence,
                            qcat_config)
        aligned_adapter, adapter_identity = ret

    if aligned_adapter is not None:
        adapter_end_position = aligned_adapter.end_query
        adapter_score = aligned_adapter.score

    return adapter_end_position, adapter_identity, adapter_score


def get_norm_socre(template, score, qcat_config):
    """
    Normalize alignment score

    :param template: Barcoding adapter
    :param score: Alignment score
    :param qcat_config: Config
    :return: Normalized score
    """
    bc_len = template.get_barcode_length(0) + template.get_barcode_length(1)
    a_len = template.get_adapter_length()
    return score * 100.0 / ((a_len - bc_len) * qcat_config.match + bc_len * qcat_config.nmatch)


def find_best_adapter_template(adapter_templates, read_sequence,
                               qcat_config):
    """
    Aligns all passed adapter templates to the read sequence returns the one
    with the highest alignment score

    :param adapter_templates: List of AdapterLayout objects
    :type adapter_templates: List
    :param read_sequence: Sequence of the read
    :type read_sequence: str
    :param qcat_config: qcatConfig object
    :type qcat_config: qcatConfig
    :return: Sequence of best adapter, alignment identity, last position of
    the aligned adapter in the read sequence,
    last position of barcode in the adapter template, length of the barcode
    :rtype: str, float, int, int, int
    """
    best_adapter_score = -1.0
    best_adapter_end_position = -1
    best_adapter_template = -1

    if not adapter_templates or not read_sequence:
        return best_adapter_template, best_adapter_end_position, best_adapter_score

    # If single adapter template is passed, convert it to list
    if not isinstance(adapter_templates, list):
        adapter_templates = [adapter_templates]

    for i, template in enumerate(adapter_templates):

        if not template.get_adapter_sequences():
            continue

        ret = eval_adapter_template(adapter_template=template,
                                    read_sequence=read_sequence,
                                    qcat_config=qcat_config,
                                    identity=False)
        adapter_end_position, _, adapter_score = ret

        adapter_score = get_norm_socre(template, adapter_score, qcat_config)

        if best_adapter_score < adapter_score:
            best_adapter_score = adapter_score
            best_adapter_template = i
            best_adapter_end_position = adapter_end_position

    return best_adapter_template, best_adapter_end_position, best_adapter_score


def build_return_dict(best_barcode, best_barcode_score,
                      best_adapter,
                      best_adapter_end,
                      exit_status,
                      trim5p=0, trim3p=0
                      ):
    """
    Builds dict containing all return values

    :param best_barcode: Barcode with the best alignment score
    :type  best_barcode: Barcode
    :param best_barcode_score: Normalized alignment score
    :type  best_barcode_score: float
    :param best_adapter: Adapter identified
    :type  best_adapter: AdapterLayout
    :param best_adapter_end: End position of adapter on the read
    :type  best_adapter_end: int
    :return: Dictionary
    """
    return_dict = {"barcode":                best_barcode,
                   "barcode_score":          best_barcode_score,
                   "adapter":                best_adapter,
                   "adapter_end":            best_adapter_end,
                   "trim5p":                 trim5p,
                   "trim3p":                 trim3p,
                   "exit_status":            exit_status
                  }

    return return_dict


def empty_return_dict():
    """
    Return empty Barcoding result dict

    :return: Empty barcode result dict
    """
    return build_return_dict(
        best_barcode=None,
        best_barcode_score=0.0,
        best_adapter=None,
        best_adapter_end=0,
        exit_status=1,
        trim3p=0,
        trim5p=0
    )


class BarcodeScanner(object):
    """
    Abstract base calss for BarcodeScanners
    """

    def __init__(self, min_quality, kit_name, kit_folder=None,
                 enable_filter_barcodes=False, scan_middle_adapter=False):
        """
        Init

        :param min_quality: Minimum barcode qulity
        :param kit_name: Name of the kit to call barcodes for
        :param kit_folder: Folder containing kit config files
        :param enable_filter_barcodes: Remove low abundance adapters in batch
        mode
        :param scan_middle_adapter: Scann full read for adapters
        """

        available_kits = adapters.populate_adapter_layouts(kit_folder)

        # Parameters
        self.min_quality = min_quality

        self.layouts = []
        self.override_kit_name = None

        self.enable_filter_barcodes = enable_filter_barcodes
        self.scan_middle_adapter = scan_middle_adapter

        # Get kets
        if kit_name and kit_name.lower() != 'auto':
            for layout in available_kits:
                if kit_name.lower() == layout.kit.lower():
                    self.layouts.append(layout)
        else:
            for layout in available_kits:
                if layout.auto_detect:
                    self.layouts.append(layout)

    @staticmethod
    def get_name():
        """
        Get name of barcode scanner algorithm

        :return: Name
        """
        raise NotImplemented("Abstract class")

    def barcode_count(self):
        """
        Number of supported barcodes

        :return: barcode number
        """
        raise NotImplemented("Abstract class")

    def scan(self, read_sequence, read_qualities, barcoding_kits,
             non_barocding_kits, qcat_config=config.qcatConfig()):
        """

        :param read_sequence:
        :param read_qualities:
        :param barcoding_kits:
        :param non_barocding_kits:
        :param qcat_config:
        :return:
        """
        raise NotImplemented("Abstract class")

    def scan_middle(self, sequence, kit_name, qcat_config):

        detected_adapters = self.get_adapters(kit_name)

        middle_barcode = self.scan(sequence[qcat_config.max_align_length:-qcat_config.max_align_length],
                                   None,
                                   detected_adapters,
                                   [],
                                   qcat_config=qcat_config)

        if not middle_barcode or middle_barcode[
            'barcode_score'] < 50.0:

            middle_barcode = self.scan(revcomp(sequence[
                                       qcat_config.max_align_length:-qcat_config.max_align_length]),
                                       None,
                                       detected_adapters,
                                       [],
                                       qcat_config=qcat_config)
            if not middle_barcode or middle_barcode[
                'barcode_score'] < 50.0:
                return False
            else:
                logging.debug(
                    "middle adapter ({}) found on rev strand with score {} at {}".format(
                        kit_name,
                        middle_barcode[
                            'barcode_score'],
                        middle_barcode[
                            'adapter_end']),
                    file=sys.stderr)
                return True
        else:
            logging.debug(
                "middle adapter ({}) found with score {} at {}".format(kit_name,
                                                                       middle_barcode[
                                                                           'barcode_score'],
                                                                       middle_barcode[
                                                                           'adapter_end']),
                file=sys.stderr)
            return True

    def detect_barcode(self,
                       read_sequence,
                       read_qualities=None,
                       qcat_config=config.qcatConfig()):

        if not self.override_kit_name:
            kits = self.layouts
        else:
            kits = self.get_adapters(self.override_kit_name)

        # Check 5' end
        align_seq_5p = extract_align_sequence(read_sequence,
                                              False,
                                              qcat_config.max_align_length)

        barcode_dict_5p = self.scan(align_seq_5p,
                                    read_qualities,
                                    kits,
                                    [],
                                    qcat_config=qcat_config)

        trim_5p = 0
        if barcode_dict_5p['adapter_end'] > 0:
            trim_5p = barcode_dict_5p['adapter_end']

        if barcode_dict_5p and barcode_dict_5p[
            'barcode_score'] < self.min_quality:
            barcode_dict_5p = empty_return_dict()

        # Check 3' end
        align_seq_5p_rc = extract_align_sequence(read_sequence,
                                                 True,
                                                 qcat_config.max_align_length)

        barcode_dict_5p_rc = self.scan(align_seq_5p_rc,
                                       read_qualities,
                                       kits,
                                       [],
                                       qcat_config=qcat_config)

        trim_3p = len(read_sequence)
        if barcode_dict_5p_rc['adapter'] and barcode_dict_5p_rc[
            'adapter_end'] > 0:
            trim_3p = trim_3p - barcode_dict_5p_rc['adapter_end']

        if barcode_dict_5p_rc and \
                barcode_dict_5p_rc['barcode_score'] < self.min_quality:
            barcode_dict_5p_rc = empty_return_dict()

        # Find best barcode
        results = [barcode_dict_5p, barcode_dict_5p_rc]

        best = None
        best_score = 0.0
        for i in range(0, len(results)):
            result = results[i]
            if result:
                if result['barcode_score'] > best_score:
                    best_score = result['barcode_score']
                    best = result

        if not best:
            best = empty_return_dict()
        else:
            min_score = 60
            if barcode_dict_5p['barcode'] and barcode_dict_5p_rc['barcode']:
                if barcode_dict_5p['barcode_score'] >= min_score and \
                        barcode_dict_5p_rc['barcode_score'] >= min_score and \
                        barcode_dict_5p['barcode'].id != barcode_dict_5p_rc['barcode'].id:
                    best = empty_return_dict()
                    best['exit_status'] = 1002

        if self.scan_middle_adapter and best['adapter'] and self.scan_middle(read_sequence, best['adapter'].kit, qcat_config):
            best = empty_return_dict()
            best['exit_status'] = 997

        best["trim5p"] = trim_5p
        best["trim3p"] = trim_3p

        if best["trim3p"] < best["trim5p"]:
            # Only happens for reads that only consist of the barcode
            best["trim5p"] = 0

        return best

    def get_adapters(self, kit_name):
        result = []
        for layout in self.layouts:
            if kit_name.lower() == layout.kit.lower():
                result.append(layout)
        return result

    def get_adapter(self, kit_name):
        for layout in self.layouts:
            if kit_name.lower() == layout.kit.lower():
                return layout

    def scan_end(self, sequence, reverse, qcat_config):
        # Check 5' end
        align_seq_5p = extract_align_sequence(sequence,
                                              reverse,
                                              qcat_config.max_align_length)

        ret = find_best_adapter_template(adapter_templates=self.layouts,
                                         read_sequence=align_seq_5p,
                                         qcat_config=qcat_config)

        best_adapter_template_index, aligned_adapter_end, best_adapter_score = ret

        return self.layouts[best_adapter_template_index], aligned_adapter_end, best_adapter_score

    def scan_ends(self, read_sequence, qcat_config):
        # 5' end
        adapter_5p, end5p, score_5p = self.scan_end(read_sequence, False, qcat_config)
        # 3' end
        adapter_3p, end3p, score_3p = self.scan_end(read_sequence, True, qcat_config)
        end3p = qcat_config.max_align_length - end3p

        if score_5p > score_3p:
            return adapter_5p, adapter_3p
        else:
            return adapter_3p, adapter_5p

    @staticmethod
    def update_kit_count(adapter, adapter_counts):
        if adapter:
            kit = adapter.kit
            if kit not in adapter_counts.keys():
                adapter_counts[kit] = 0
            adapter_counts[kit] += 1
        else:
            if "none" not in adapter_counts.keys():
                adapter_counts["none"] = 0
            adapter_counts["none"] += 1

    @staticmethod
    def get_most_abundant_kits(adapter_counts):
        if not adapter_counts:
            return None
        return sorted(adapter_counts.items(), key=operator.itemgetter(1), reverse=True)[0][0]

    def detect_kit(self, read_sequences, qcat_config):

        adapter_counts = {}

        adapter_ends = []

        for read_sequence in read_sequences:
            adapter_1, adapter_2 = self.scan_ends(read_sequence, qcat_config)
            self.update_kit_count(adapter_1, adapter_counts)
            # self.update_kit_count(adapter_2, adapter_counts)
            # self.update_kit_count(adapter, adapter_counts)

            # adapter_ends.append((end5p, end3p))

        kit_name = self.get_most_abundant_kits(adapter_counts)

        return kit_name, adapter_ends

    @staticmethod
    def update_barcode_count(result, barcode_count):
        if result and result['barcode']:
            if result['barcode'].id not in barcode_count:
                barcode_count[result['barcode'].id] = 0
            barcode_count[result['barcode'].id] += 1
        else:
            if "0" not in barcode_count:
                barcode_count["0"] = 0
            barcode_count["0"] += 1

    @staticmethod
    def get_valid(barcode_counts, min_perc=0.20):
        max_barcode_count = 0
        for n in barcode_counts.values():
            max_barcode_count = max(max_barcode_count, n)

        min_count = int(max_barcode_count * min_perc)

        valid_barcodes = []
        for bc, count in barcode_counts.items():
            if count > min_count:
                valid_barcodes.append(bc)

        return valid_barcodes

    def filter_barcodes(self, barcode_count, results):
        valid_barcodes = self.get_valid(barcode_count, 0.05)

        for i, barcode in enumerate(results):
            if barcode and barcode['barcode'] and barcode['barcode'].id not in valid_barcodes:
                results[i] = empty_return_dict()
        return results

    def detect_barcode_batch(self, read_sequences, read_qualities=[None],
                             qcat_config=config.qcatConfig()):
        # barcode_count = [0] * 1000
        barcode_count = {}

        kit_name, _ = self.detect_kit(read_sequences, qcat_config)
        results = []

        self.override_kit_name = kit_name
        for read_sequence, read_quality in zip(read_sequences, read_qualities):
            result = self.detect_barcode(read_sequence, read_quality, qcat_config)
            self.update_barcode_count(result, barcode_count)
            results.append(result)

        self.override_kit_name = None

        if self.enable_filter_barcodes:
            results = self.filter_barcodes(barcode_count, results)

        return results
