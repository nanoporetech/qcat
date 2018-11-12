from __future__ import print_function

import logging
import time
import sys
import os

from six.moves import zip
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
    parser.add_argument("-t", "--tsv", dest="TSV", default=None, required=False)
    parser.add_argument("-n", "--name", dest="NAME", default="-", required=False)
    parser.add_argument("-d", "--dataset", dest="DATASET", default="-", required=False)
    parser.add_argument("-s", "--summary", action="store_true", dest="SUMMARY", required=False)
    parser.add_argument("-g", "--guppy_summary", dest="GUPPY", default=None, required=False)
    parser.add_argument("-i", "--get-incorrect", dest="INCORRECT", default=None, required=False)
    args = parser.parse_args(argv)

    return args


def _parse_reads_info(comment):
    """
    Parses ONT FASTQ header information

    :param comment: FASTQ comment
    :return: dict
    """
    single_read_info = {}
    if comment:
        for pair in comment.split(" "):
            cols = pair.split('=')
            if len(cols) == 2:
                single_read_info[cols[0]] = cols[1]
    return single_read_info


def from_guppy(fastq):
    with open(fastq) as fh:

        for line in fh:
            cols = line.split("\t")
            if len(cols) > 3:
                yield cols[1].split(','), cols[2]


def get_col_number(line, pattern):
    cols = line.split('\t')
    for i, name in enumerate(cols):
        if name.strip().lower() == pattern.strip().lower():
            return i
    return None


def from_tsv(fastq):
    if fastq == '-':
        fastq = "/dev/stdin"
    with open(fastq) as fh:
        for line in fh:
            if line.startswith('name'):
                comment_col = get_col_number(line, 'comment')
                name_col = get_col_number(line, 'name')
                barcode_col = get_col_number(line, 'barcode')
                score_col = get_col_number(line, 'score')
                continue
            cols = line.strip().split("\t")
            if len(cols) > 3:
                info = _parse_reads_info(cols[comment_col])
                truebc = info['truebc'].split(",")
                yield cols[name_col], truebc, cols[barcode_col], float(cols[score_col])


def from_stdin(fastq, guppy_calls=None):
    for names, comments, sequences, qualites in cli.iter_fastx(fastq, True, 1):

        name, comment, sequence, quality = next(zip(names, comments, sequences, qualites))

        info = _parse_reads_info(comment)
        truebc = info['truebc'].split(",")
        bc = "none"
        if guppy_calls and name in guppy_calls:
            bc, score = guppy_calls[name]
            if bc == "unclassified":
                bc = "none"
        else:
            if info['barcode']:
                bc = info['barcode']
            score = -1

        yield name, truebc, bc, score


def main(argv=sys.argv[1:]):
    """
    Basic command line interface to qcat.

    :param argv: Command line arguments
    :type argv: list
    :return: None
    :rtype: NoneType
    """
    args = parse_args(argv=argv)

    correct = 0
    incorrect = 0
    fn = 0
    fp = 0

    summary = args.SUMMARY
    if not summary:
        print("program", "dataset", "read_id", "truebc", "bc", "result", "score", sep='\t')

    incorrect_only = args.INCORRECT

    guppy_calls = {}
    if args.GUPPY:
        with open(args.GUPPY) as fh:
            fh.readline()
            for line in fh:
                cols = line.split("\t")
                guppy_calls[cols[0]] = (cols[1].replace("barcode0", "").replace("barcode", ""), cols[5])

    if args.TSV:
        generator = from_tsv(args.TSV)
    else:
        generator = from_stdin(args.FASTQ, guppy_calls)

    for id, truebc, bc, score in generator:
        if bc in truebc:
            correct += 1
            if not summary:
                print(args.NAME, args.DATASET, id, ",".join(truebc), bc, "CORRECT", score, sep='\t')
        else:
            if bc == "none":
                if truebc != ["none"]:
                    fn += 1
                    if not summary:
                        print(args.NAME, args.DATASET, id, ",".join(truebc), bc, "FN", score, sep='\t')
            else:
                if truebc == ["none"]:
                    fp += 1
                    if not summary:
                        print(args.NAME, args.DATASET, id, ",".join(truebc), bc, "FP", score, sep='\t')
                else:
                    incorrect += 1
                    if not summary:
                        print(args.NAME, args.DATASET, id, ",".join(truebc), bc, "INCORRECT", score, sep='\t')


    total = correct + incorrect + fn + fp
    if total == 0:
        print("No reads found", file=sys.stderr)
        sys.exit(1)
    correct_perc = correct * 100.0 / total
    incorrect_perc = (incorrect + fp) * 100.0 / total
    fn_perc = fn * 100.0 / total

    if summary:
        print(args.NAME, args.DATASET, total, correct, incorrect, fn, fp, correct_perc,
                  incorrect_perc, fn_perc, sep="\t")


if __name__ == '__main__':

    main()
