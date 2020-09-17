
# idemuxCPP - inline barcode demultiplexing
idemuxCPP is a command line tool designed to demultiplex paired-end fastq files from
[QuantSeq-Pool](https://www.lexogen.com/quantseq-pool-sample-barcoded-3mrna-sequencing/).

idemuxCPP can demultiplex based on i7, i5 and i1 inline barcodes. While this tool
can generally be used to demultiplex on any barcodes (as long as they are correctly supplied
and in the fastq header), it best performs when used in combination with
[Lexogen indices](https://www.lexogen.com/indexing/12nt-dual-indexing-kits/), as it
will correct common sequencing errors in the sequenced barcodes. This will allow you
to retain more reads from your sequencing experiment, while minimizing cross contamination.


idemuxCPP use is permitted under the following [licence](LICENCE).

idemuxCPP is a direct translation of the python tool idemux (https://github.com/lexogen-tools/idemux)
in order to decrease the runtime. It is 2 times faster than the python version.

**General usage:**
```
    idemuxCPP [-h] --r1 READ1 --r2 READ2 [--sample-sheet SAMPLE_SHEET] --out OUTPUT_DIR
           [--i1-start I1_START] [--i5-rc] [-v]
```

**Run idemuxCPP:**
```
    idemuxCPP --r1 read_1.fastq.gz --r2 read_2.fastq.gz --sample-sheet samples.csv --out /some/output/path --i1-start pos_in_read_2
```

## Features
* FASTQ file demultiplexing based on i7, i5 or i1 barcodes
* Correction of barcode sequencing errors to maximize read yield (only works
  with [Lexogen 12 nt UDIs](https://www.lexogen.com/indexing/12nt-dual-indexing-kits/),
  that have been sequenced at least 8 nt.
* Reverse complementation in case the i5 index has been sequenced as reverse complement


## Getting started
To get stated with demultiplexing you need to:

1. [Install idemuxCPP](#1-installation)
2. [Prepare a sample sheet csv](#2-preparing-the-sample-sheet)
3. [Run idemuxCPP](#3-running-idemuxcpp)

## 1. Installation
dependencies:

* compiler supporting C++11 standard and OpenMP
* [boost C++ library](http://www.boost.org/) version >= 1.55.0 (install the development versions of the following libraries (or install all e.g. in Ubuntu via package `libboost-all-dev`)

  * libboost_filesystem
  * libboost_system
  * libboost-test (only required if you want to compile unit tests)
* [zlib](https://zlib.net/) (e.g. zlib1g-dev in Ubuntu)
* [gengetopt](https://www.gnu.org/software/gengetopt/gengetopt.html)


**Windows 10 64bit binary**
For windows you do not need to install any dependencies (they are included in the package).
Simply download the pre-compiled windows binary from here [windows binary](win_solution/iDemux_win10_64bit.zip).
Extract the zip file. To execute the tool press `windows+r`, enter `cmd`, `cd C:\\location_of_the_extracted_zip_file\bin` and execute `.\\idemuxCPP`


**From Source**

To configure, compile and install execute the following commands on your command line:
```
    ./configure [--help for additional configuration options]
    make
    make install
```

**From Linux Package**
<table><thead><tr>
<th> Debian </th>
<th> Ubuntu </th>
<th> Fedora </th>
</tr></thead><tbody><tr>
<td style="vertical-align:top">
<details><summary>Debian_10</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Debian_10/i386/idemuxcpp_0.1.2-1_i386.deb"> idemuxcpp - 0.1.2 - 32 bit</a></p>
<p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Debian_10/amd64/idemuxcpp_0.1.2-1_amd64.deb"> idemuxcpp - 0.1.2 - 64 bit</a></p>
</details>
<details><summary>Debian_9.0</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Debian_9.0/i386/idemuxcpp_0.1.2-1_i386.deb"> idemuxcpp - 0.1.2 - 32 bit</a></p>
<p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Debian_9.0/amd64/idemuxcpp_0.1.2-1_amd64.deb"> idemuxcpp - 0.1.2 - 64 bit</a></p>
</details></td>
<td style="vertical-align:top">
<details><summary>xUbuntu_20.04</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/xUbuntu_20.04/amd64/idemuxcpp_0.1.2-1_amd64.deb"> idemuxcpp - 0.1.2 - 64 bit</a></p></details>
<details><summary>xUbuntu_19.04</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/xUbuntu_19.04/amd64/idemuxcpp_0.1.2-1_amd64.deb"> idemuxcpp - 0.1.2 - 64 bit</a></p></details>
</details><details><summary>xUbuntu_18.04</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/xUbuntu_18.04/amd64/idemuxcpp_0.1.2-1_amd64.deb"> idemuxcpp - 0.1.2 - 64 bit</a></p>
</details><details><summary>xUbuntu_16.04</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/xUbuntu_16.04/i386/idemuxcpp_0.1.2-1_i386.deb"> idemuxcpp - 0.1.2 - 32 bit</a></p>
<p><a href="https://download.opensuse.org/repositories/home:/Lexogen/xUbuntu_16.04/amd64/idemuxcpp_0.1.2-1_amd64.deb"> idemuxcpp - 0.1.2 - 64 bit</a></p></details>
</td>
<td style="vertical-align:top">
<details><summary>Fedora 32</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Fedora_32/x86_64/idemuxcpp-0.1.2-147.1.x86_64.rpm"> idemuxcpp - 0.1.2 - 64 bit</a></p></details>
</td>
</tr></tbody></table>



idemuxCPP will also soon be available via bioconda!


## 2. Preparing the sample sheet
In order to run idemuxCPP on your QuantSeq-Pool data you first need to prepare a 
[csv file](https://en.wikipedia.org/wiki/Comma-separated_values).
We call this csv a sample sheet and it specifies which barcodes correspond to each
sample.

This is a necessity as the software needs to know into which bins reads should be
sorted during demultiplexing. A sample sheet can easily be generated by filling in an
excel spreadsheet and exporting it as csv.


Example sample sheet (i7, i5 and i1 demuliplexing):
```
    sample_name,i7,i5,i1
    sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
    sample_1,AAAATCCCAGTT,CCCCTAAACGTT,AAAATCCCAGTT
    sample_2,GAAAATTTACGC,GCCCCTTTCAGA,GAAAATTTACGC
    sample_3,AAACTAACTGTC,CCCATCCATGTA,AAACTAACTGTC
```

A sample sheet consists of 4 columns and  always starts with the header illustrated
above. 'Sample_name' values will be used as output file names, while the
sequences specified in i7,i5 & i1 will be used for demultiplexing.

Therefore, only specific, unique unambiguous combinations of sample names and barcodes are
allowed. This means using duplicated or ambiguous combinations will result in an error.
However, idemuxCPP will do its best to tell you where the problem lies, once this happens.


**In brief the rules are:**

1. Sample names need to be unique.
2. Barcode combinations need to be unique.
3. i7 and/or i5 indices have to be used consistently within the csv file.
   i7 and/or i5 indices need to be either present for all samples or none at all.
4. In contrast to i7/i5 indices, i1 indices can be used for a subset of samples in the csv file.
5. Absence of a barcode needs to be indicated by an empty field (no value between
   comas ``,,``).
6. If your i5 has been sequenced as reverse complement, *do not* enter the reverse
   complement sequences in the sample sheet. Use the ``--i5-rc`` option!


See [below](#sample-sheet-examples). for more showcases of sample/barcode combinations that are *allowed* or
*disallowed*.


## 3. Running idemuxCPP
Once you have installed the tool you can run it by typing ``idemuxCPP`` in the terminal.

idemuxCPP accepts the following arguments:
```
    required arguments:
      --r1 READ1                   path to gzipped read 1 fastq file
      --r2 READ2                   path to gzipped read 2 fastq file
      --sample-sheet CSV           csv file describing sample names, and barcode combinations
      --out OUTPUT_DIR             where to write the output files

    optional arguments:
      --i5-rc                      when the i5 barcode has been sequenced as reverse complement.
                                   make sure to always use non-reverse complement sequences in the sample sheet
      --i1_start POS               start position of the i1 index (1-based) on read 2 (default: 11)
      -v, --version                show program's version number and exit
      -h, --help                   show help message and exit
```

Example commands:
```
    # demultiplexes read 1 and 2 into the folder 'demux'
    idemuxCPP --r1 read_1.fastq.gz --r2 read_2.fastq.gz --sample-sheet samples.csv --out demux

    # demultiplexing assuming the i1 barcode starts at the first base
    idemuxCPP --r1 read_1.fastq.gz --r2 read_2.fastq.gz --sample-sheet samples.csv --out demux --i1_start 1

    # demultiplexing assuming i5 is present as reverse complement in the fastq header
    # if he i5 has been sequenced as reverse complement use this option and provide
    # the NON reverse complement sequences in the sample sheet.
    idemuxCPP --r1 read_1.fastq.gz --r2 read_2.fastq.gz --sample-sheet samples.csv --out demux
```

After a successful completed run idemuxCPP will write summary report to the output folder
('demultipexing_stats.tsv').

## Technicalities
When you run idemuxCPP the following will happen:

* It will check if your sample sheet is okay. See [here](#sample-sheet-examples) for examples

* It will check the fastq header for barcodes and expects them in the following format:
    ```
    single index (i7 or i5): @NB502007:379:HM7H2BGXF:1:11101:24585:1069 1:N:0:TCAGGTAANNTT

    dual index (i7 and i5): @NB502007:379:HM7H2BGXF:1:11101:24585:1069 1:N:0:TCAGGTAANNTT+NANGGNNCNNNN
    ```
* Reads that cannot be demultiplexed will be written to undetermined_R{1/2}.fastq.gz

* When you demultiplex based on i1 inline barcodes, the a successful recognized barcode
  sequence will be cut out and removed from read 2. This is a design choice and will leave
  you with the 10 nt UMI + the nucleotides that potentially follow the i1 barcode
  (or don't).

This allows you to:

1. Use other software, such as UMI_tools to deal with the 10nt UMI if desired
2. To demuliplex lanes where QuantSeq-Pool has been pooled with other libraries and read
   2 has been sequenced longer than the actual barcode.

If you sequenced i5 as a reverse complement, make sure to not fill in reverse complement
barcodes into the sample sheet, but to use the ``--i5-rc`` parameter.

## Help
If you are demuliplexing a large number of samples (more than 500) you might encounter the
following error:

* ``OSError: [Errno 24] Too many open files``

This error occurs because most OS have a limit on how many files can be opened and
written to at the ame time. In order to temporarily increase the limit run:
```
    # multiply your sample number*2 (as data is paired end)
    # then round to the next multiple of 1024
    ulimit -n the_number_above
```
If you are looking for a permanent solution you can change your ulimit values
[this way](https://access.redhat.com/solutions/61334).

In case you experience any issues with this software please open an issue describing your
problem. Make sure to post the version of the tool you are running (``-v, --version``)
and your os.

## Sample sheet examples
*This is allowed:*
```
    # demultiplexing via full i7, i5, i1
    sample_name,i7,i5,i1
    sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
    sample_1,AAAATCCCAGTT,CCCCTAAACGTT,AAAATCCCAGTT

    # demultiplexing via full i7, i5 and sparse i1
    sample_name,i7,i5,i1
    sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
    sample_1,AAAATCCCAGTT,CCCCTAAACGTT,

    # demultiplexing via full i7, i5
    sample_name,i7,i5,i1
    sample_0,AAAACATGCGTT,CCCCACTGAGTT,
    sample_1,AAAATCCCAGTT,CCCCTAAACGTT,

    # demultiplexing via full i7, no i5 and sparse i1
    sample_name,i7,i5,i1
    sample_0,AAAACATGCGTT,,AAAACATGCGTT
    sample_1,AAAATCCCAGTT,,

    # demultiplexing via full i7 only
    sample_name,i7,i5,i1
    sample_0,AAAACATGCGTT,,
    sample_1,AAAATCCCAGTT,,

    # demultiplexing via full i5 and i1
    sample_name,i7,i5,i1
    sample_0,,CCCCACTGAGTT,AAAACATGCGTT
    sample_1,,CCCCTAAACGTT,AAAATCCCAGTT

    # demultiplexing via full i5 and sparse i1
    sample_name,i7,i5,i1
    sample_0,,CCCCACTGAGTT,AAAACATGCGTT
    sample_1,,CCCCTAAACGTT,

    # demultiplexing via full i5
    sample_name,i7,i5,i1
    sample_0,,CCCCACTGAGTT,
    sample_1,,CCCCTAAACGTT,

    # demultiplexing via full i1
    sample_name,i7,i5,i1
    sample_0,,,AAAACATGCGTT
    sample_1,,,AAAATCCCAGTT
```

*This is not allowed:*
```
    # missing i1 column (or any other)
    sample_name,i7,i5,
    sample_0,AAAACATGCGTT,CCCCACTGAGTT
    sample_1,AAAATCCCAGTT,CCCCTAAACGTT

    # duplicated barcode combination
    sample_name,i7,i5,i1
    sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
    sample_1,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT

    # duplicated sample names
    sample_name,i7,i5,i1
    sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
    sample_0,AAAATCCCAGTT,CCCCTAAACGTT,AAAATCCCAGTT

    # mixed, potentially ambiguous indexing (full i7 and sparse i5, i1)
    sample_name,i7,i5,i1
    sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
    sample_1,AAAATCCCAGTT,,AAAATCCCAGTT
    sample_2,GAAAATTTACGC,GCCCCTTTCAGA,GAAAATTTACGC
    sample_3,AAACTAACTGTC,,AAACTAACTGTC

    # mixed, potentially ambiguous indexing indexing (no i7, sparse i5 & i1)
    sample_name,i7,i5,i1
    sample_0,,CCCCACTGAGTT,
    sample_1,,,AAAATCCCAGTT

    # mixed, potentially ambiguous indexing indexing (sparse i7, full i5 & i1)
    sample_name,i7,i5,i1
    sample_0,,CCCCACTGAGTT,AAAACATGCGTT
    sample_1,AAAATCCCAGTT,CCCCTAAACGTT,AAAATCCCAGTT
    sample_2,,GCCCCTTTCAGA,GAAAATTTACGC
    sample_3,AAACTAACTGTC,CCCATCCATGTA,AAACTAACTGTC

    # missing comma separator
    sample_name,i7,i5,i1
    sample_0,AAAACATGCGTTCCCCACTGAGTT,AAAACATGCGTT

    # no barcodes
    sample_name,i7,i5,i1
    sample_0,,,

    # wrong column headers
    wrong_col_name,i7,i5,i1
    sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
```

&copy; Lexogen GmbH, 2020

