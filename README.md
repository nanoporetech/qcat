![.](ONT_logo.png  "Oxford Nanopore Technologies")
******************

# qcat

`qcat` is Python command-line tool for demultiplexing Oxford Nanopore reads from FASTQ files. It accepts basecalled FASTQ files and splits the reads into into separate FASTQ files based on their barcode. Qcat makes the demultiplexing algorithms used in albacore/guppy and EPI2ME available to be used locally with FASTQ files. Currently qcat implements the EPI2ME algorithm. In the next version we will add the albacore/guppy algorithm.

|Algorithm | Status |Description  |
|--|--|--|
|Guppy/Albacore | not yet available | Using this mode, qcat will produce demultiplexing results identical to Albacore. |
| EPI2ME | available | Using this mode, qcat will produce demultiplexing results identical to EPI2ME's demultiplexing workflow |

In addition, qcat supports demultiplexing of datasets prepared using dual or combinatorial barcodes (see "How to run qcat?" section for more information). 

If you want to demultiplex your reads during basecalling, please use [albacore](https://github.com/nanoporetech/albacore).

Change log
------------
**v1.1.0**
- Add RAB214

**v1.0.7**
- Fix PBC and PBK adapter trimming

**v1.0.6**
- Fixed issue with file names in dual barcoding mode
- Resolved problem when using –tsv and -b at the same time
- Min. read length filter added
- Fixed YAMLLoadWarning

Requirements
------------
* Linux or MacOS

Quick start
-----------
For the vast majority of datasets, it will be sufficient to run qcat using default parameters:
```bash
$ qcat -f <fastq_file> -b <output folder>
```
After qcat finished, please check to summary output to verify that a barcode was assigned to most of the reads.

Installation 
------------
**Conda (recommended)** 

[![Anaconda-Server Badge](https://anaconda.org/bioconda/qcat/badges/version.svg)](https://anaconda.org/bioconda/qcat)

To install qcat using conda, make sure you have [Miniconda3](https://conda.io/miniconda.html) installed and the [bioconda](https://bioconda.github.io/#install-conda) channels set up. After bioconda is set up, you can install qcat as follows:
```bash
$ conda install qcat
```

**PIP** 

[![PyPI version](https://badge.fury.io/py/qcat.svg)](https://badge.fury.io/py/qcat)

If you want to install qcat using pip, please make sure to use python3.
```bash
$ pip install --user qcat
```
To install qcat using pip for all users on a system, run (requires root permissions):
```bash
$ pip install qcat
```

**Docker**

If you have [docker](https://docs.docker.com/install/) available on your computer, you can run qcat using the following command without any prior installation:
```bash
$ docker run -ti -v `pwd`:`pwd` -w `pwd` quay.io/biocontainers/qcat:1.0.0--py_0 qcat -f ./input_file.fastq -b ./output_folder/
```
When running qcat using the above command, the input and output folders have to be in the current working directory. If you want to read from or write to locations other than your working directory, please adjust the `-v` parameter accordingly.

**Manually install from source**

[![GitHub version](https://badge.fury.io/gh/nanoporetech%2Fqcat.svg)](https://badge.fury.io/gh/nanoporetech%2Fqcat)

To install qcat manually, please make sure you have python3 and git available, and install as follows:
```bash
$ git clone https://github.com/nanoporetech/qcat.git
$ cd qcat
$ python3 setup.py install
```

How to run qcat?
----------------

**Demultiplexing multiple FASTQ files from a folder**
```bash
$ cat input_folder/*.fastq | qcat -b <output folder>
```
**Demultiplexing a single FASTQ file**
```bash
$ qcat -f <fastq_file> -b <output folder>
```
**Demultiplexing dual/combinatorial barcoding datasets**
```bash
$ qcat --dual -f <fastq_file> -b <output folder>
```

What does the output look like?
--------------------------------

### Demultiplexing
Qcat will give you a summary of what barcodes where found and how many on the command line. For example:
```bash
Adapters detected in 191 of 193 reads
         LWB001    191: |  ################### |  98.96 %
           none      2: |                      |   1.04 %
Barcodes detected in 191 of 193 adapters
      barcode01    191: |  ################### |  98.96 %
           none      2: |                      |   1.04 %
Demultiplexing finished in 0.69s
```
In addition, you will get an output folder with a single FASTQ file per barcode:
```bash
$ ls -R tmp/output/
barcode01.fastq none.fastq
```
These files can be used for downstream analysis.


Frequently asked questions
--------------------------

### Why are there different demultiplexing modes?
Currently Albacore and EPI2ME use different algorithms for demultiplexing. To provide maximum compatibility, qcat will support both algorithms. Eventually, both algorithms will be consolidated.  

### What is the difference between porechop and qcat?
Qcat works similarity to porechop. However, porechop was recently deprecated.  

### When I start qcat without any command-line parameters, it starts but does not terminate.
When invoked without `-f/--fastq`, qcat will wait for input from standard in. Press Ctrl+C to stop qcat and either rerun wiht `-f/--fastq`:
```bash
$ qcat -f input_file.fastq -b output_folder/
```
or make sure to use a pipe to provide the input FASTQ file. E.g.: 
```bash
$ cat input_file.fastq | qcat -b output_folder/
```
### Qcat expects a single input file, but I got several FASTQ files from basecalling.

A simple option is to pipe all FASTQ files into qcat:
```bash
$ cat basecalled/*.fastq | qcat -b output_folder/
```
Alternatively, you can concatenate all FASTQ files into a single file first:
```bash
$ cat basecalled/*.fastq > single_file.fastq
```
and use `-f/--fastq`:
```bash
$ qcat -f single_file.fastq -b output_folder/
```
### I want to use Albacore's algorithm for demultiplexing, but get a warning saying "Demultiplexing mode guppy currently not supported. Falling back to epi2me."

Currently, Albacore's demultiplexing algorithm is only supported when running the qcat docker image. In the next version, we will support other ways of running Guppy/Albacore demultiplexing as well.

### Can qcat trim adapters without barcodes?

Currently qcat only supports demultiplexing. However, in the future we plan to add a trimming only mode as well.

### What is "simple" mode?

In simple mode qcat will only search for barcode sequences without trying to determine which barcoding kit was used. ONLY USE THIS MODE FOR TESTING/DEBUGGING PURPOSE. Albacore and EPI2ME mode will be faster and give you better results.

### What barcoding kits does qcat support?

| Kit  | Description  |
|--|--|
| Auto  | Auto detect kit |
| RBK001  | Rapid barcoding kit |
| RBK004 | Rapid barcoding kit v4 |
| NBD103/NBD104 | Native barcoding kit with barcodes 1-12 |
| NBD114  | Native barcoding kit with barcodes 13-24 |
| NBD104/NBD114 | Native barcoding kit with barcodes 1-24 |
| PBC001 | PCR Barcoding Kit with 12 barcodes |
| PBC096 | PCR Barcoding Kit with 96 barcodes |
| RPB004/RLB001  | Rapid PCR Barcoding Kit (SQK-RPB004) and Rapid Low Input by PCR Barcoding Kit |
| PBK004/LWB001 | Low Input by PCR Barcoding Kit |
| RAB204 | 16S Rapid Amplicon Barcoding Kit with 12 Barcodes  |
| VMK001 | Voltrax Barcoding Kit with 4 barcodes |


Full usage
------------
```
usage: qcat [-h] [-V] [-l LOG] [--quiet] [-f FASTQ] [-b BARCODE_DIR]
            [-o OUTPUT] [--min-score MIN_QUAL] [--detect-middle] [-t THREADS]
            [--tsv] [--trim]
            [-k {Auto,PBC096,RBK004,NBD104/NBD114,PBK004/LWB001,RBK001,RAB204,VMK001,PBC001,NBD114,NBD103/NBD104,DUAL,RPB004/RLB001}]
            [--list-kits] [--guppy | --epi2me | --dual | --simple]
            [--no-batch] [--filter-barcodes]
            [--simple-barcodes {standard,extended}]

Python command-line tool for demultiplexing Oxford Nanopore reads from FASTQ files

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
  -l LOG, --log LOG     Print debug information
  --quiet               Don't print summary

General settings:
  -f FASTQ, --fastq FASTQ
                        Barcoded read file
  -b BARCODE_DIR, --barcode_dir BARCODE_DIR
                        If specified, qcat will demultiplex reads to this
                        folder
  -o OUTPUT, --output OUTPUT
                        Output file trimmed reads will be written to (default:
                        stdout).
  --min-score MIN_QUAL  Minimum barcode score. Barcode calls with a lower
                        score will be discarded. Must be between 0 and 100.
                        (default: 60)
  --detect-middle       Search for adapters in the whole read
  -t THREADS, --threads THREADS
                        Number of threads. Only works with in guppy mode
  --tsv                 Prints a tsv file containing barcode information each
                        read to stdout.
  --trim                Remove adapter and barcode sequences from reads.
  -k {Auto,PBC096,RBK004,NBD104/NBD114,PBK004/LWB001,RBK001,RAB204,VMK001,PBC001,NBD114,NBD103/NBD104,DUAL,RPB004/RLB001}, --kit {Auto,PBC096,RBK004,NBD104/NBD114,PBK004/LWB001,RBK001,RAB204,VMK001,PBC001,NBD114,NBD103/NBD104,DUAL,RPB004/RLB001}
                        Sequencing kit. Specifying the correct kit will
                        improve sensitivity and specificity and runtime
                        (default: auto)
  --list-kits           List all supported kits

Demultiplexing modes:
  --guppy               Use Guppy's demultiplexing algorithm (default: false)
  --epi2me              Use EPI2ME's demultiplexing algorithm (default: true)
  --dual                Use dual barcoding algorithm
  --simple              Use simple demultiplexing algorithm. Only looks for
                        barcodes, not for adapter sequences. Use only for
                        testing purposes!

EPI2ME options (only valid with --epi2me):
  --no-batch            Don't use information from multiple reads for kit
                        detection (default: false)
  --filter-barcodes     Filter rare barcode calls when run in batch mode

Simple options (only valid with --simple):
  --simple-barcodes {standard,extended}
                        Use 12 (standard) or 96 (extended) barcodes for
                        demultiplexing
```

**Licence and Copyright**

© 2018 Oxford Nanopore Technologies Ltd.

`qcat` is distributed under the terms of the Mozilla Public License 2.0.

**Research Release**

Research releases are provided as technology demonstrators to provide early
access to features or stimulate Community development of tools. Support for
this software will be minimal and is only provided directly by the developers.
Feature requests, improvements, and discussions are welcome and can be
implemented by forking and pull requests. However much as we would
like to rectify every issue and piece of feedback users may have, the 
developers may have limited resource for support of this software. Research
releases may be unstable and subject to rapid iteration by Oxford Nanopore
Technologies.