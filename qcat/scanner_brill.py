import os

import pkg_resources

import numpy as np

from qcat import config
from qcat.scanner import BarcodeScanner
from qcat.scanner_base import empty_return_dict, build_return_dict

try:
    from brill.model import load_model
    from brill import seq_util
    from brill import encoding
    use_brill = True
except ImportError as e:
    use_brill = False


class BarcodeScannerBrill(BarcodeScanner):

    def __init__(self, min_quality=None, kit_folder=None, kit=None, enable_filter_barcodes=False, scan_middle_adapter=False):

        if min_quality is None:
            min_quality = 10

        super(BarcodeScannerBrill, self).__init__(min_quality,
                                                  kit,
                                                  kit_folder=kit_folder,
                                                  enable_filter_barcodes=enable_filter_barcodes,
                                                  scan_middle_adapter=scan_middle_adapter)

        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'

        # Load all available brill models
        self.brill_models = {}
        for layout in self.layouts:
            ret = self.load_model(layout)
            if ret:
                self.brill_models[layout.kit] = ret

    def load_model(self, adapter):
        brill_model = self.find_model(adapter.model)
        if not brill_model:
            return None

        # Load model:
        return load_model(brill_model), adapter.model_len

    @staticmethod
    def get_name():
        return "brill"

    @staticmethod
    def run_brill(model, model_len, read_sequence, end5p=None, end3p=None):
        if not end5p:
            end5p = model_len
        if not end3p:
            end3p = model_len

        x = encoding.build_predictor_matrix(read_sequence[end5p - model_len:end5p],
                                            read_sequence[-end3p:][:model_len])

        prediction = model.predict_proba(x, verbose=0)

        return BarcodeScannerBrill.process_prediction(prediction)

    @staticmethod
    def find_model(path):
        if not path:
            return None
        if os.path.exists(path):
            return path
        else:
            resource_model = pkg_resources.resource_filename(__name__, os.path.join("resources/models", path))
            if os.path.exists(resource_model):
                return resource_model
            else:
                raise RuntimeError("Couldn't locate model {}".format(path))

    def detect_barcode(self,
                       read_sequence,
                       read_qualities=None,
                       qcat_config=config.qcatConfig()):

        result = empty_return_dict()

        if not self.override_kit_name:
            detected_adapter, _ = self.scan_ends(read_sequence, qcat_config)
            override_kit_name = detected_adapter.kit
        else:
            override_kit_name = self.override_kit_name

        if override_kit_name in self.brill_models:
            brill_model, brill_model_len = self.brill_models[override_kit_name]

            if len(read_sequence) > brill_model_len:
                bc, qscore = self.run_brill(brill_model, brill_model_len, read_sequence)

                barcode = None
                adapter = None
                if bc and bc != 'none':
                    adapter = self.get_adapter(override_kit_name)
                    barcode = adapter.get_barcode_set()[int(bc) - 1]

                if bc and qscore >= self.min_quality:
                    result = build_return_dict(
                        best_barcode=barcode,
                        best_barcode_score=qscore,
                        best_adapter=adapter,
                        best_adapter_end=90,
                        trim5p=90,
                        trim3p=90,
                        exit_status=0
                    )

        return result

    @staticmethod
    def process_prediction(p):
        """ Process CNN prediction, return barcode and q-score.

        :param p: CNN prediction.
        :param min_qual: Minimum quality for classified barcodes.
        :return: barcode and qscore
        :rtype: tuple

        """
        # Select winning barcode:
        bc = np.argmax(p[0])
        # Calculate error probability and q-score:
        prob_error = 1 - p[0][bc]
        qscore = seq_util.prob_to_phred(prob_error, max_q=200)

        if bc == 0:
            # Unclassified bin predicted:
            bc = None
        return bc, qscore


