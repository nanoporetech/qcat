try:
    import ConfigParser
except ImportError:
    import configparser as ConfigParser
import parasail


class qcatConfig:

    def __init__(self, config_path=None):

        self._match = 5
        self._nmatch = -1
        self._mismatch = -2
        self._gap_open = 2
        self._gap_extend = 2

        self._max_align_length = 150

        self._extracted_barcode_extension = 11
        self._barcode_context_length = 11

        self._matrix = None
        self.update_matrix()

        self._matrix_barcode = parasail.matrix_create("ATGCN", 1, -1)

        if config_path is not None:
            self.read(config_path)

    @property
    def matrix_barcode(self):
        """
        Scoring matrix used for aligning barcodes

        :return: parasail scoring matrix
        """
        return self._matrix_barcode


    @property
    def match(self):
        """
        Match score for adapter/barcode alignments

        :return: float
        """
        return self._match

    @match.setter
    def match(self, value):
        """
        Match score for adapter/barcode alignments

        :param value: float
        :return: None
        """
        self._match = abs(value)
        self.update_matrix()

    @property
    def nmatch(self):
        """
        Match score for adapter/barcode alignments

        :return: float
        """
        return self._nmatch

    @nmatch.setter
    def nmatch(self, value):
        """
        Match score for adapter/barcode alignments

        :param value: float
        :return: None
        """
        self._nmatch = abs(value)
        self.update_matrix()

    @property
    def mismatch(self):
        """
        Mismatch score for adapter/barcode alignments

        :return: float
        """
        return self._mismatch

    @mismatch.setter
    def mismatch(self, value):
        """
        Mismatch score for adapter/barcode alignments

        :param value: float
        :return: None
        """
        self._mismatch = -1 * abs(value)
        self.update_matrix()

    @property
    def gap_open(self):
        """
        Gap open for adapter/barcode alignments

        :return: float
        """
        return self._gap_open

    @gap_open.setter
    def gap_open(self, value):
        """
        Gap open for adapter/barcode alignments

        :param value: float
        :return: None
        """
        self._gap_open = abs(value)

    @property
    def gap_extend(self):
        """
        Gap extend for adapter/barcode alignments

        :return: float
        """
        return self._gap_extend

    @gap_extend.setter
    def gap_extend(self, value):
        """
        Gap extend for adapter/barcode alignments

        :param value: float
        :return: None
        """
        self._gap_extend = abs(value)

    @property
    def max_align_length(self):
        """
        Number of bp of the read that will be scanned for an adapter

        :return: int
        """
        return self._max_align_length

    @max_align_length.setter
    def max_align_length(self, value):
        """
        Number of bp of the read that will be scanned for an adapter

        :param value: int
        :return: None
        """
        self._max_align_length = value

    @property
    def extracted_barcode_extension(self):
        """
        Number of bp the detected barcode region from the read will be
        extended before aligning barcodes

        :return: int
        """
        return self._extracted_barcode_extension

    @extracted_barcode_extension.setter
    def extracted_barcode_extension(self, value):
        """
        Number of bp the detected barcode region from the read will be
        extended before aligning barcodes

        :param value: int
        :return: None
        """
        self._extracted_barcode_extension = value

    @property
    def barcode_context_length(self):
        """
        Number of bp upstream and downstream of the barcoded t
        hat will be added to the barcode sequence
        when aligning to the read.

        :return: int
        """
        return self._barcode_context_length

    @barcode_context_length.setter
    def barcode_context_length(self, value):
        """
        Number of bp upstream and downstream of the barcoded t
        hat will be added to the barcode sequence
        when aligning to the read.

        :param value: int
        :return: None
        """
        self._barcode_context_length = value

    def write(self, out_config_path):
        """
        Write to ini file

        :param out_config_path: ini file
        :return: None
        """
        default_config = ConfigParser.RawConfigParser()
        default_config.add_section('qcat')
        default_config.set('qcat',
                           'gap_open',
                           self.gap_open)
        default_config.set('qcat',
                           'gap_extend',
                           self.gap_extend)
        default_config.set('qcat',
                           'match',
                           self.match)
        default_config.set('qcat',
                           'mismatch',
                           self.mismatch)
        default_config.set('qcat',
                           'max_align_length',
                           self.max_align_length)
        default_config.set('qcat',
                           'extracted_barcode_extension',
                           self.extracted_barcode_extension)
        default_config.set('qcat',
                           'barcode_context_length',
                           self.barcode_context_length)

        with open(out_config_path, 'wb') as configfile:
            default_config.write(configfile)

    def update_matrix(self):
        """
        Create new parasail scoring matrix. 'N' is used as wildcard character
        for barcodes and has its own match parameter (0 per default).
        'X' is used as wildcard character for modified bp as in the 16S
        sequencing adapter.

        :return: None
        """
        self.matrix = parasail.matrix_create("ATGCNX", self.match, self.mismatch)

        pointers = [4, 11, 18, 25, 28, 29, 30, 31, 32]
        for i in pointers:
            self.matrix.pointer[0].matrix[i] = self.nmatch

        pointers = [5, 12, 19, 26, 33, 35, 36, 37, 38, 39, 40]
        for i in pointers:
            self.matrix.pointer[0].matrix[i] = 0

    def read(self, config_path):
        """
        Read from ini file

        :param config_path: ini path
        :type config_path: str
        :return: None
        """
        config = ConfigParser.RawConfigParser()
        config.read(config_path)

        self.gap_open = config.getint('qcat', 'gap_open')
        self.gap_extend = config.getint('qcat', 'gap_extend')
        self.match = config.getint('qcat', 'match')
        self.mismatch = config.getint('qcat', 'mismatch')

        self.max_align_length = config.getint('qcat',
                                              'max_align_length')
        self.extracted_barcode_extension = config.getint('qcat',
                                                         'extracted_barcode_extension')
        self.barcode_context_length = config.getint('qcat',
                                                    'barcode_context_length')


def get_default_config():
    """
    Return default config object which can be modified and passed to
    detect_barcode_* functions

    :return: Config object
    :rtype CuecCatConfig
    """
    return qcatConfig()
