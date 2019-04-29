from __future__ import print_function

import logging
import time
import sys
import os
import six

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from argparse import ArgumentParser, RawDescriptionHelpFormatter, ArgumentTypeError

from qcat import __version__, adapters, config
from qcat import scanner
from qcat.adapters import Barcode
from qcat.scanner import get_modes, factory, get_kits_info, get_kits


def get_mode(args):
    if args.MODE_GUPPY:
        logging.warning("Guppy/albacore barcoding is not yet supported. "
                        "Falling back to epi2me.")
        return "epi2me"
    if args.MODE_EPI2ME:
        return "epi2me"
    if args.MODE_SIMPLE:
        return "simple"
    if args.MODE_DUAL:
        return "dual"

    return "epi2me"


def check_minqual_arg(x):
    x = float(x)
    if x < 0.0 or x > 100.0:
        raise ArgumentTypeError("Minimum quality must be a value between 0 and 100.")
    return x


def check_kit_arg(x):
    x = str(x)
    if x.lower() == "dual":
        raise ArgumentTypeError(
            "-k dual and --kit dual are not supported any more. "
            "Please use --dual instead.")
    return x


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
    defaults = config.get_default_config()
    usage = "Python command-line tool for demultiplexing Oxford Nanopore reads from FASTQ files"
    parser = ArgumentParser(description=usage,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-V", '--version',
                        action='version',
                        version='%(prog)s ' + __version__)
    parser.add_argument('-l', "--log",
                        dest="log",
                        type=str,
                        default="INFO",
                        help="Print debug information")
    parser.add_argument("--quiet",
                        dest="QUIET",
                        action='store_true',
                        help="Don't print summary")

    general_group = parser.add_argument_group('General settings')
    general_group.add_argument("-f", "--fastq",
                        type=str,
                        dest="fastq",
                        help="Barcoded read file")
    general_group.add_argument('-b', "--barcode_dir",
                        dest="barcode_dir",
                        type=str,
                        default=None,
                        help="If specified, qcat will demultiplex reads "
                             "to this folder")
    general_group.add_argument('-o', "--output",
                        dest="output",
                        type=str,
                        default=None,
                        help="Output file trimmed reads will be written to "
                             "(default: stdout).")
    general_group.add_argument("--min-score",
                        dest="min_qual",
                        type=check_minqual_arg,
                        default=None,
                        help="Minimum barcode score. Barcode calls with a lower "
                             "score will be discarded. "
                             "Must be between 0 and 100. (default: 60)")
    general_group.add_argument("--detect-middle",
                        dest="DETECT_MIDDLE",
                        action='store_true',
                        help="Search for adapters in the whole read")
    general_group.add_argument('-t', "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads. Only works with in guppy mode")
    general_group.add_argument("--min-read-length",
                               dest="min_length",
                               type=int,
                               default=100,
                               help="Reads short than <min-read-length> after "
                                    "trimming will be discarded.")
    general_group.add_argument("--tsv",
                        dest="tsv",
                        action='store_true',
                        help="Prints a tsv file containing barcode information "
                             "each read to stdout.")
    general_group.add_argument("--trim",
                        dest="TRIM",
                        action='store_true',
                        help="Remove adapter and barcode sequences from reads.")
    general_group.add_argument("-k", "--kit",
                               dest="kit",
                               choices=get_kits(),
                               type=check_kit_arg,
                               default="auto",
                               help="Sequencing kit. Specifying the correct kit "
                                    "will improve sensitivity and specificity and "
                                    "runtime (default: auto)")
    general_group.add_argument("--list-kits",
                               dest="list_kits",
                               action='store_true',
                               help="List all supported kits")

    # Mode
    mode_group = parser.add_argument_group("Demultiplexing modes")
    group = mode_group.add_mutually_exclusive_group()

    group.add_argument("--guppy",
                        dest="MODE_GUPPY",
                        action = 'store_true',
                        help="Use Guppy's demultiplexing algorithm (default: false)")
    group.add_argument("--epi2me",
                       dest="MODE_EPI2ME",
                       action='store_true',
                       help="Use EPI2ME's demultiplexing algorithm (default: true)")
    group.add_argument("--dual",
                       dest="MODE_DUAL",
                       action='store_true',
                       help="Use dual barcoding algorithm")
    group.add_argument("--simple",
                       dest="MODE_SIMPLE",
                       action='store_true',
                       help="Use simple demultiplexing algorithm. Only looks for "
                            "barcodes, not for adapter sequences. Use only for testing "
                            "purposes!")

    # EPI2ME
    epi2me_group = parser.add_argument_group('EPI2ME options (only valid with --epi2me)')
    epi2me_group.add_argument("--no-batch",
                              dest="nobatch",
                              action='store_true',
                              help="Don't use information from multiple reads for kit detection (default: false)")
    epi2me_group.add_argument("--filter-barcodes",
                              dest="FILTER_BARCODES",
                              action='store_true',
                              help="Filter rare barcode calls when run in batch mode")

    simple_group = parser.add_argument_group('Simple options (only valid with --simple)')
    simple_group.add_argument("--simple-barcodes",
                              dest="SIMPLE_BARCODES",
                              # choices=['standard', 'extended'],
                              default="standard",
                              help="Use 12 (standard) or 96 (extended) barcodes for demultiplexing")

    # parser.add_argument("--adapter-sequences",
    #                     type=str,
    #                     dest="adapter_yaml",
    #                     default=None,
    #                     help="YAML file or folder containing YAML files that"
    #                          "describe the Adapter sequences.")

    # parser.add_argument('--no-header',
    #                     action='store_true',
    #                     dest='no_header',
    #                     help="Do not write header to output tsv file."
    #                          "Useful if results from different runs should be"
    #                          "combined")

    args = parser.parse_args(argv)

    return args


def extract_fastx_comment(header):
    """
    Split FASTA/Q header into name and comment

    :param header: FASTA/Q header line
    :return: name, comment
    :rtype: tuple
    """
    cols = header.replace("\t", " ").split(" ")
    name = cols[0]
    comment = None
    if len(cols) > 1:
        comment = " ".join(cols[1:])

    return name, comment


def is_fastq(filename):
    if not filename:
        logging.debug("Reading from stdin. Assuming FASTQ format. "
                     "Please use -f/--fasta to run with FASTA files.")
        return True
    c = None
    with open(filename) as f:
        c = f.read(1)

    if c == "@":
        return True
    elif c == ">":
        return False
    else:
        raise ValueError("Invalid input file. File must start with "
                           "'@' or '>'. Current file starts with: " + c)


def iter_fastx(reads_fx, fastq, batchsize):
    """
    Return iterator for FASTA/Q file

    :param reads_fx: filename of FASTX file
    :return: None
    """
    # batch = []

    names = []
    comments = []
    seqs = []
    quals = []

    # fastq = is_fastq(reads_fx)
    if fastq:

        if reads_fx:
            logging.debug("Reading from file {}".format(reads_fx))
            f = open(reads_fx)
        else:
            logging.debug("Reading from stdin")
            f = sys.stdin

        try:
            for title, seq, qual in FastqGeneralIterator(f):
                name, comment = extract_fastx_comment(title)

                # batch.append((name, comment, seq, qual))
                names.append(name)
                comments.append(comment)
                seqs.append(seq)
                quals.append(qual)

                if len(names) >= batchsize:
                    yield names, comments, seqs, quals
                    # batch = []
                    names = []
                    comments = []
                    seqs = []
                    quals = []
        except ValueError as e:
            logging.error(e.message)
            sys.exit(1)

        if reads_fx:
            logging.debug("Closing file {}".format(reads_fx))

                # yield name, comment, seq, qual

    else:
        with open(reads_fx) as f:
            for title, seq in SimpleFastaParser(f):
                name, comment = extract_fastx_comment(title)

                # batch.append((name, comment, seq, None))
                names.append(name)
                comments.append(comment)
                seqs.append(seq)
                quals.append(None)

                if len(names) >= batchsize:
                    yield names, comments, seqs, quals
                    # batch = []
                    names = []
                    comments = []
                    seqs = []
                    quals = []

    if len(names) > 0:
        # yield batch
        yield names, comments, seqs, quals


def get_output_file(output_files, out_folder, barcode_dict, fastq):

    bc_id = "none"
    if barcode_dict['barcode']:
        bc_id = barcode_dict['barcode'].name
        if bc_id:
            bc_id = bc_id.replace('/', '_')

    if bc_id not in output_files.keys():
        if fastq:
            output_files[bc_id] = open(
                os.path.join(out_folder, bc_id + ".fastq"), "w")
        else:
            output_files[bc_id] = open(
                os.path.join(out_folder, bc_id + ".fasta"), "w")
    return output_files[bc_id]


def write_to_file(trimmed_output_file, output_files, out_folder,
                  name, comment, sequence, quality, fastq, result):

    if not comment:
        comment = ""
    if not quality:
        quality = ""

    if out_folder:
        out_file = get_output_file(output_files, out_folder, result, fastq)
        if fastq:
            print("@" + name + " " + comment, sequence, "+", quality, sep="\n",
                  file=out_file)
        else:
            print(">" + name + " " + comment, sequence, sep="\n",
                  file=out_file)

    # if trimmed_output_file:
    else:

        bc = "none"
        if result["barcode"]:
            bc = str(result['barcode'].id)

        comment = "{} barcode={}".format(comment, bc)

        if fastq:
            print("@" + name + " " + comment, sequence, "+", quality, sep="\n",
                  file=trimmed_output_file)
        else:
            print(">" + name + " " + comment, sequence, sep="\n",
                  file=trimmed_output_file)


def close_files(output_files):
    for f in output_files.values():
        f.close()


def adapter_found(adapter_dist, adapter):
    adapterid = "none"
    if adapter:
        adapterid = adapter.kit

    if adapterid not in adapter_dist.keys():
        adapter_dist[adapterid] = 0
    adapter_dist[adapterid] += 1


def barcode_found(barcode_dist, barcode):
    bcid = "none"
    if barcode:
        bcid = barcode.name

    if bcid not in barcode_dist.keys():
        barcode_dist[bcid] = 0
    barcode_dist[bcid] += 1


def print_barcode_hist(barcode_dist, adapter_dist, total_reads):
    adapters_detected = 0
    for key, value in six.iteritems(adapter_dist):
        if key != "none":
            adapters_detected += value
    logging.info("Adapters detected in %d of %d reads" % (adapters_detected, total_reads))
    for adapterid in sorted(adapter_dist):
        perc = adapter_dist[adapterid] * 100.0 / total_reads
        logging.info("%15s %6d: | %20s | %6s %%" %(adapterid, adapter_dist[adapterid],
                                     int(perc / 5) * "#", "{:.2f}".format(perc)))

    barcodes_detected = 0
    for key, value in six.iteritems(barcode_dist):
        if key != "none":
            barcodes_detected += value
    logging.info("Barcodes detected in %d of %d adapters" % (barcodes_detected, total_reads))
    for bcid in sorted(barcode_dist):
        perc = barcode_dist[bcid] * 100.0 / total_reads
        logging.info("%15s %6d: | %20s | %6s %%" % (bcid, barcode_dist[bcid],
                                    int(perc / 5) * "#", "{:.2f}".format(perc)))


def write_multiplexing_result(barcode_dict, comment, name, sequence, tsv):
    """
    Writes results to TSV (--tsv) file.

    :param barcode_dict: detected barcode
    :param comment: FASTQ comment
    :param name:  FASTQ name
    :param sequence: FASTQ sequence
    :param tsv: TSV path
    :return: None
    """
    if barcode_dict['barcode']:
        kit_name = None
        if barcode_dict['adapter']:
            kit_name = barcode_dict['adapter'].kit

        if tsv:
            print(name,
                  len(sequence),
                  barcode_dict['barcode'].id,
                  barcode_dict['barcode_score'],
                  kit_name,
                  barcode_dict['adapter_end'],
                  comment,
                  sep="\t")
    else:
        if tsv:
            print(name,
                  len(sequence),
                  "none",
                  "-1",
                  "none",
                  "-1",
                  comment,
                  sep='\t')


def qcat_cli(reads_fq, kit, mode, nobatch, out,
               min_qual, tsv, output, threads, trim, adapter_yaml, quiet, filter_barcodes, middle_adapter, min_read_length,
               qcat_config):
    """
    Runs barcode detection for each read in the fastq file
    and print the read name + the barcode to a tsv file

    :param reads_fq: Path to fastq file
    :type reads_fq: str
    :param no_header: Print header or not for output tsv file
    :type no_header: bool
    :param kit: Which kit was used for sequencing
    :type kit: str
    :param simple: Triggers "simple" detection mode where only
    :type simple: bool
    barcodes will be aligned to the read. Adapter detection is skipped.
    :param nobatch: Switch of "batch" mode (First the kit and the barcodes used
    in the sample are detected from a small subset of read. Then, only these
    adapter/barcodes will be used for barcode detection)
    :type nobatch: bool
    :param qcat_config: qcatConfig object
    :type qcat_config: qcatConfig
    :return: None
    """

    total_reads = 0
    skipped_reads = 0
    barcode_dist = {}
    adapter_dist = {}

    notrimming = not trim

    detector = factory(mode=mode,
                       kit=kit,
                       min_quality=min_qual,
                       kit_folder=adapter_yaml,
                       enable_filter_barcodes=filter_barcodes,
                       scan_middle_adapter=middle_adapter,
                       threads=threads)

    if tsv:
        print("name", "length", "barcode", "score", "kit",
              "adapter_end", "comment", sep="\t")

    fastq = is_fastq(reads_fq)
    output_files = {}

    trimmed_output_file = sys.stdout
    if output:
        trimmed_output_file = open(output, "w")

    if out:
        if not os.path.exists(out):
            os.makedirs(out)

    batch_size = 4000
    if nobatch:
        batch_size = 1

    for names, comments, seqs, quals in iter_fastx(reads_fq, fastq, batch_size):
        # Detect adapter/barcode
        if nobatch:
            results = [detector.detect_barcode(read_sequence=seqs[0],
                                               read_qualities=quals[0],
                                               qcat_config=qcat_config)]
        else:
            results = detector.detect_barcode_batch(read_sequences=seqs,
                                                    read_qualities=quals,
                                                    qcat_config=qcat_config)
        for name, comment, sequence, quality, result in zip(names,
                                                            comments,
                                                            seqs,
                                                            quals,
                                                            results):
            total_reads += 1

            if not notrimming:
                trim_5p = result["trim5p"]
                trim_3p = result["trim3p"]
                sequence = sequence[trim_5p:trim_3p]
                if quality:
                    quality = quality[trim_5p:trim_3p]

            if len(sequence) < min_read_length:
                skipped_reads += 1
                continue

            # Record which adapter/barcode was found
            barcode_found(barcode_dist, result['barcode'])
            adapter_found(adapter_dist, result['adapter'])

            # Write tsv result file
            write_multiplexing_result(result,
                                      comment,
                                      name,
                                      sequence,
                                      tsv)
            # Write FASTQ/A files
            if out or not tsv:
                write_to_file(trimmed_output_file,
                              output_files,
                              out,
                              name,
                              comment,
                              sequence,
                              quality,
                              fastq,
                              result)

    if out:
        close_files(output_files)

    if not quiet:
        print_barcode_hist(barcode_dist, adapter_dist, total_reads)
        if skipped_reads > 0:
            logging.info("{} reads were skipped due to the min. length filter.".format(skipped_reads))

    if output:
        trimmed_output_file.close()


def barcodes_from_fasta(filename):
    barcodes = []
    for names, comments, seqs, quals in iter_fastx(filename, False, 1):
        for name, comment, sequence in zip(names, comments, seqs):
            barcodes.append(Barcode("BC{}".format(name),
                                    int(name),
                                    sequence,
                                    True))

    return barcodes


def main(argv=sys.argv[1:]):
    """
    Basic command line interface to qcat.

    :param argv: Command line arguments
    :type argv: list
    :return: None
    :rtype: NoneType
    """
    try:
        args = parse_args(argv=argv)

        qcat_config = config.get_default_config()

        numeric_level = getattr(logging, args.log.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: %s' % args.log.upper())
        logging.basicConfig(level=numeric_level, format='%(message)s')

        if args.list_kits:
            kits = get_kits_info()
            for kit in sorted(kits.keys()):
                if kit != "auto" and kit != "DUAL":
                    logging.info("{:<30}{}".format(kit, kits[kit]))
            return

        mode = get_mode(args)
        kit = args.kit
        if mode == "simple":
            kit = args.SIMPLE_BARCODES

        start = time.time()
        qcat_cli(reads_fq=args.fastq,
                 # no_header=args.no_header,
                 kit=kit,
                 mode=mode,
                 nobatch=args.nobatch,
                 out=args.barcode_dir,
                 min_qual=args.min_qual,
                 tsv=args.tsv,
                 output=args.output,
                 threads=args.threads,
                 trim=args.TRIM,
                 adapter_yaml=None,#=args.adapter_yaml, #args.barcode_fa,
                 quiet=args.QUIET,
                 filter_barcodes=args.FILTER_BARCODES,
                 middle_adapter=args.DETECT_MIDDLE,
                 min_read_length=args.min_length,
                 qcat_config=qcat_config)
        end = time.time()

        if not args.QUIET:
            logging.info("Demultiplexing finished in {0:.2f}s".format(end - start))
    except IOError as e:
        logging.error(e)
    except ValueError as e:
        logging.error(e)


if __name__ == '__main__':

    main()
