from __future__ import print_function

import logging

try:
    import pysam
except ImportError as e:
    logging.warning(
        "Could not load pysam. To use qcat-eval-full please install pysam.")
import math
import sys
try:
    import mappy as mp
except ImportError as e:
    logging.warning(
        "Could not load pysam. To use qcat-eval-full please install mappy.")

import os
import glob

from argparse import ArgumentParser, RawDescriptionHelpFormatter


from qcat import config
from qcat.scanner import BarcodeScanner

LOOKUP = []


def qstring_to_phred(quality):
    """ Compute standard phred scores from a quality string. """
    qscores = []
    if quality is not None:
        qscores = [ord(q) - 33 for q in quality]
    return qscores


for q in range(100):
    LOOKUP.append(pow(10, -.1 * q))


def _compute_mean_qscore(scores):
    """ Returns the phred score corresponding to the mean of the probabilities
    associated with the phred scores provided.

    :param scores: Iterable of phred scores.

    :returns: Phred score corresponding to the average error rate, as
        estimated from the input phred scores.
    """
    if not scores:
        return 0.0
    sum_prob = 0.0
    for val in scores:
        sum_prob += LOOKUP[val]
    mean_prob = sum_prob / len(scores)
    return -10.0 * math.log10(mean_prob)



def parse_args(argv):
    """
    Commandline parser

    :param argv: Command line arguments
    :type argv: List
    :param defaults: qcat config object that holds default values passed to
    ArgumentParser
    :type defaults: qcatConfig
    :return: None
    """
    usage = "Command line interface to barcode detection library"
    parser = ArgumentParser(description=usage,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument(dest="FASTQ", nargs="?", default=None)
    parser.add_argument("-t", "--tsv", dest="TSV", required=True)
    parser.add_argument("-n", "--name", dest="NAME", required=True)
    parser.add_argument("-g", "--genomes", dest="GENOMES", required=True)
    parser.add_argument("--max", dest="MAX", default=-1, required=False)
    # parser.add_argument("-s", "--summary", action="store_true", dest="SUMMARY", required=False)
    # parser.add_argument("-g", "--guppy_summary", dest="GUPPY", default=None, required=False)
    args = parser.parse_args(argv)

    return args


def _parse_reads_info(comment):
    """
    Parse ONT FASTQ header

    :param comment: str
    :return: dict
    """
    single_read_info = {'runid': 'unknown'}
    if comment:
        for pair in comment.split(" "):
            cols = pair.split('=')
            if len(cols) == 2:
                single_read_info[cols[0]] = cols[1]
    return single_read_info


def create_reference(mapping, reference_file, fasta_path):
    ref_count = 0
    with open(reference_file, "w") as ref:
        for g in glob.glob(fasta_path):
            with pysam.FastxFile(g) as fh:
                id = os.path.splitext(os.path.basename(g))[0]
                for entry in fh:
                    if id in mapping:
                        ref_count += 1
                        name = ','.join([str(x) for x in mapping[id]])

                        print(">" + str(mapping[id][0]) + " " + entry.name,
                              entry.comment, file=ref)
                        print(entry.sequence, file=ref)
    if ref_count != len(mapping.keys()):
        raise RuntimeError(
            "Couldn't find all references, please check mappings!")
    return ref_count


def get_truebc(a, sequence, min_aln_perc=0.85, adapter_len=200):
    #             print(bc, score)
    hits = list(a.map(sequence))

    truebc = None
    truebc_assignment = None
    if len(hits) == 0:
        truebc_assignment = "Not_found"
    else:
        if len(hits) == 1:
            hit = hits[0]
            if hit.blen > min_aln_perc * (len(sequence) - adapter_len):
                truebc = hit.ctg
                truebc_assignment = "Unique"
            else:
                truebc = hit.ctg
                truebc_assignment = "Partial"
        else:
            #             for hit in hits:
            #                 print(hit)

            truebcs = set([hit.ctg for hit in hits])
            truebc = ",".join(truebcs)
            if len(truebcs) == 1:
                truebc_assignment = "Fragmented"
            else:
                truebc = "Multiple"
                truebc_assignment = "Multiple"
    #             print(truebc, truebc_assignment)

    #             print(truebc, truebc_assignment)
    #             print("---")

    return truebc, truebc_assignment


def main(argv=sys.argv[1:]):
    """
    Basic command line interface to qcat.

    :param argv: Command line arguments
    :type argv: list
    :return: None
    :rtype: NoneType
    """
    args = parse_args(argv=argv)

    reference_file = args.GENOMES

    a = mp.Aligner(reference_file, preset="map-ont")  # load or build index

    if not a: raise Exception("ERROR: failed to load/build index")

    output_file = open(args.TSV, "w")

    detector = BarcodeScanner.factory(kit="NBD10X", min_quality=0)

    print("dataset", "id", "length", "qscore", "score", "barcode", "truebc",
          "truebc_assignment", "middle_adapter", "non_matching_bc", "sequence",
          sep='\t', file=output_file)
    name = args.NAME
    f = args.FASTQ
    i = 0
    print(name, f)
    with pysam.FastxFile(f) as fh:
        for entry in fh:
            i += 1

            info = _parse_reads_info(entry.comment)
            if 'truebc' in info:
                truebc = info['truebc']
                truebc_assignment = "Pre_assigned"
            else:
                truebc, truebc_assignment = get_truebc(a, entry.sequence)

            c = config.get_default_config()
            result = detector.detect_barcode(entry.sequence,
                                             qcat_config=c)

            middle_adapter = detector.scan_middle(entry.sequence, "NBD10X",
                                                  c)

            q_score = _compute_mean_qscore(qstring_to_phred(entry.quality))
            bc = "none"
            score = 0.0
            if result['barcode']:
                bc = result['barcode'].id
                score = result['barcode_score']

            no_seq = False
            seq = "-"
            if not no_seq:
                seq = entry.sequence

            print(name, entry.name, len(entry.sequence), q_score, score, bc,
                  truebc, truebc_assignment, middle_adapter,
                  result['exit_status'] == 1002, seq,
                  sep='\t', file=output_file)

            #             if i == 10000:
            #                 break
            if i % 1000 == 0:
                print(".", end="")
            if i % 50000 == 0:
                print("")

            if i == int(args.MAX):
                print("Reach maximumg, stopping.")
                break
    print("")

    output_file.close()

    print("finished")


if __name__ == '__main__':

    main()
