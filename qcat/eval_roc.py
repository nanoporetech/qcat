from __future__ import print_function

import logging

try:
    import pandas as pd
except ImportError as e:
    logging.warning("Could not load pandas. To use qcat-roc please install pandas.")
import sys

from argparse import ArgumentParser, RawDescriptionHelpFormatter

from qcat import cli


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
    # parser.add_argument("-t", "--tsv", dest="TSV", default=None, required=False)
    # parser.add_argument("-n", "--name", dest="NAME", default="-", required=False)
    # parser.add_argument("-d", "--dataset", dest="DATASET", default="-", required=False)
    # parser.add_argument("-s", "--summary", action="store_true", dest="SUMMARY", required=False)
    # parser.add_argument("-g", "--guppy_summary", dest="GUPPY", default=None, required=False)
    args = parser.parse_args(argv)

    return args


def summary(truebc_col, bc_col, score_col, min_score):
    counts = {'total': 0, 'correct': 0, 'incorrect': 0, 'true_negative': 0,
              'false_negative': 0, 'false_positive': 0, 'score': min_score}
    for truebc, bc, score in zip(truebc_col, bc_col, score_col):
        if score < min_score:
            bc = "none"

        truebc = [str(truebc)]
        if bc == "none":
            if truebc == ["none"] or truebc == "none":
                counts["true_negative"] += 1
            else:
                counts["false_negative"] += 1
        else:
            if truebc == ["none"] or truebc == "none":
                counts["false_positive"] += 1
            else:
                if str(bc) in truebc:
                    counts['correct'] += 1
                else:
                    counts['incorrect'] += 1
        counts['total'] += 1

    return counts


def create_roc(input_file):
    data = pd.read_csv(input_file, sep="\t")

    col_names = ['total', 'correct', 'incorrect', 'true_negative',
                 'false_negative', 'false_positive', 'score', 'program',
                 'dataset']
    my_df = pd.DataFrame(columns=col_names)

    for q in range(0, 1001):
        counts = summary(data['truebc'], data['bc'], data['score'], q / 10.0)

        counts['program'] = data['program'][0]
        counts['dataset'] = data['dataset'][0]
        my_df = my_df.append(counts, ignore_index=True)
    return my_df


def main(argv=sys.argv[1:]):
    """
    Basic command line interface to qcat.

    :param argv: Command line arguments
    :type argv: list
    :return: None
    :rtype: NoneType
    """
    args = parse_args(argv=argv)

    if args.FASTQ:
        df = create_roc(args.FASTQ)
    else:
        df = create_roc(sys.stdin)

    df.to_csv("/dev/stdout", index=False, sep='\t')
    # print(df.to_string(index=False))

if __name__ == '__main__':

    main()
