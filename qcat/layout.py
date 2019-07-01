import re

from collections import namedtuple

BarcodePosition = namedtuple("BarcodePosition", "start end length")


class AdapterLayout:
    """
    Represents the layout of adapter. Usually one kit has a specific
    adapter layout
    """

    def __init__(self, kit, sequence, barcode_set_1, barcode_set_2, description, auto_detect=False, model=None, model_len=None, name=None, trim_offset=0):
        """
        Init

        :param name: Name of adapter/kit
        :param layout_id: Unique ID
        :param sequence: Descriptive name for the adapter
        :param barcode_set_1: Set of barcodes supported by the kit
        :param barcode_set_2: Set of secondary barcodes supported (double barconding)
        :param kits: List of kits using this adapter
        :param description: Free text description
        """
        self.kit = kit
        self.auto_detect = auto_detect
        self.model = model
        self.model_len = model_len
        self.trim_offset = trim_offset
        if not name:
            self.name = self.kit
        else:
            self.name = name

        self.sequence = sequence.upper()
        if not self.sequence or re.findall("[^ATGCNX]", self.sequence):
            raise RuntimeError("Invalid adapter sequence: {}".
                               format(self.sequence))

        self.barcode_set_1 = barcode_set_1
        self.barcode_set_2 = barcode_set_2
        self.barcode_count = 0

        for barcode_set in [self.barcode_set_1, self.barcode_set_2]:
            if barcode_set:
                self.barcode_count += 1

        self.barcode_pos_1 = BarcodePosition(-1, -1, 0)
        self.barcode_pos_2 = BarcodePosition(-1, -1, 0)

        if self.barcode_set_1:
            self.barcode_pos_1 = self.get_placeholder_pos(self.sequence, 0)
            for barcode in self.barcode_set_1:
                if len(barcode.sequence) != self.barcode_pos_1.length:
                    raise RuntimeError("Adapter length does not match "
                                       "place holder length: {}, {}".
                                       format(len(barcode.sequence),
                                              self.barcode_pos_1.length))

        if self.barcode_set_2:
            self.barcode_pos_2 = self.get_placeholder_pos(self.sequence, 1)
            for barcode in self.barcode_set_2:
                if len(barcode.sequence) != self.barcode_pos_2.length:
                    raise RuntimeError("Adapter length does not match "
                                       "place holder length: {}, {}".
                                       format(len(barcode.sequence),
                                              self.barcode_pos_2.length))

        self.description = description

    @staticmethod
    def get_placeholder_pos(adapter_template, index=0):
        """
        Returns start and end position of the placehoder
        region for the barcode in the adapter

        :param adapter_template: adapter sequence with Ns
        instead of barcode sequence
        :return: Start and end position of barcode placeholder
        in adapter template
        :rtype: int, int
        """
        start = -1
        end = -1
        length = 0

        matches = list(re.finditer('N+', adapter_template))

        if matches:
            if len(matches) > index:
                start = matches[index].start(0)
                end = matches[index].end(0) - 1
                length = end - start + 1

        return BarcodePosition(start, end, length)

    def __repr__(self):
        return {"Kit": self.kit,
                "Description": self.description}.__repr__()

    def get_barcode_end(self, index=0):
        """
        Get the position of the last barcode bp in the adapter

        :param index: 0 for single barcoding, 0 or 1 for double barcoding
        :return:
        """
        if index == 0:
            return self.barcode_pos_1.end
        elif index == 1:
            return self.barcode_pos_2.end
        else:
            raise RuntimeError("Invalid barcode index: {}. Must be 0 or 1 "
                               "(for double barcoding)".format(index))

    def get_barcode_length(self, index=0):
        """
        Length of barcode in template. Usually the same for all kits

        :param index: 0 for single barcoding, 0 or 1 for double barcoding
        :return:
        """
        if index == 0:
            return self.barcode_pos_1.length
        elif index == 1:
            return self.barcode_pos_2.length
        else:
            raise RuntimeError("Invalid barcode index: {}. Must be 0 or 1 "
                               "(for double barcoding)".format(index))

    def get_adapter_sequences(self, barcode_seq=None):
        """
        Returns the sequence of the adapter. If barcode_seq is omitted,
        the region containing the barcode will be masked with Ns.

        :param barcode_seq: sequence of the barcode
        :return: str
        """
        if barcode_seq:
            seq_5p = self.sequence[:self.barcode_pos_1.start]
            seq_3p = self.sequence[self.barcode_pos_1.end + 1:]
            return seq_5p + barcode_seq + seq_3p
        else:
            return self.sequence

    def get_full_adapter_sequences(self, context=None):
        """
        Returns all possible adapter sequences including barcodes

        :param context:
        :return: Barcode, Adapter sequence
        :rtype: Barcode, List
        """
        if self.barcode_count == 0:
            yield None, self.sequence
        elif self.barcode_count == 1:
            seq_5p = self.sequence[:self.barcode_pos_1.start]
            seq_3p = self.sequence[self.barcode_pos_1.end + 1:]
            for barcode in self.barcode_set_1:
                full_sequence = seq_5p + barcode.sequence + seq_3p
                if context:
                    full_sequence = seq_5p[
                                    :-context] + barcode.sequence + seq_3p[
                                                                    context:]
                yield barcode, full_sequence

    def get_adapter_length(self):
        """
        Length of full adapter sequence including barcode

        :return: int
        """
        return len(self.get_adapter_sequences())

    def get_barcode_set(self, index=0):
        """
        Get a list of barcodes supported by the kit

        :param index: 0 for single barcoding, 0 or 1 for double barcoding
        :return:
        """
        if index == 0:
            return self.barcode_set_1
        elif index == 1:
            return self.barcode_set_2
        else:
            raise RuntimeError("Invalid barcode index: {}. Must be 0 or 1 "
                               "(for double barcoding)".format(index))

    def get_upstream_context(self, n, index=0):
        """
        Return n bp upstream of the barcode

        :param index: 0 for single barcoding, 0 or 1 for double barcoding
        :type n: int
        :param n: Number of bp
        :type n: int
        :return:
        """
        if index == 0:
            barcode = self.barcode_pos_1
        elif index == 1:
            barcode = self.barcode_pos_2
        else:
            raise RuntimeError("Invalid barcode index: {}. Must be 0 or 1 "
                               "(for double barcoding)".format(index))

        if barcode.end > -1:
            upstream_start = max(0, barcode.start - n)
            return self.sequence[upstream_start:barcode.start]
        else:
            return ""

    def get_downstream_context(self, n, index=0):
        """
        Return n bp downstream of the barcode

        :param index: 0 for single barcoding, 0 or 1 for double barcoding
        :type n: int
        :param n: Number of bp
        :type n: int
        :return:
        """
        if index == 0:
            barcode = self.barcode_pos_1
        elif index == 1:
            barcode = self.barcode_pos_2
        else:
            raise RuntimeError("Invalid barcode index: {}. Must be 0 or 1 "
                               "(for double barcoding)".format(index))

        if barcode.end > -1:
            downstream_end = min(self.get_adapter_length(),
                                 barcode.end + n + 1)
            return self.sequence[barcode.end + 1:downstream_end]
        else:
            return ""

    def is_doulbe_barcode(self):
        """
        Returns true if the kit has two barcodes per adapter
        end of the read

        :return: if kit has two barcodes per adapter
        :rtype: bool
        """
        return self.barcode_set_2 is not None
