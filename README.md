
# idemuxCPP - inline barcode demultiplexing
[![GitHub release](https://img.shields.io/github/release/Lexogen-Tools/idemuxcpp.svg)](https://github.com/Lexogen-Tools/idemuxcpp/releases)
[![Build Status](https://app.travis-ci.com/Lexogen-Tools/idemuxcpp.svg?branch=master)](https://app.travis-ci.com/github/Lexogen-Tools/idemuxcpp)
[![Conda](https://img.shields.io/conda/v/bioconda/idemuxcpp.svg)](https://anaconda.org/bioconda/idemuxcpp)
[![Conda Downloads](https://anaconda.org/bioconda/idemuxcpp/badges/downloads.svg)](https://anaconda.org/bioconda/idemuxcpp)

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


## Getting started
To get stated with demultiplexing you need to:

1. [Install idemuxCPP](#1-installation)
2. [Prepare a sample sheet csv](#2-preparing-the-sample-sheet)
3. [Run idemuxCPP](#3-running-idemuxcpp)

## 1. Installation
dependencies:

* compiler supporting C++11 standard and OpenMP
* [boost C++ library](http://www.boost.org/) version >= 1.55.0 (install the development versions of the following libraries (or install all e.g. in Ubuntu via package `libboost-all-dev`)

  * libboost-filesystem
  * libboost-system
  * libboost-iostreams
  * libboost-test (only required if you want to compile unit tests)
* [zlib](https://zlib.net/) (e.g. zlib1g-dev in Ubuntu)
* [gengetopt](https://www.gnu.org/software/gengetopt/gengetopt.html)


**Windows 10 64bit binary**
For windows you do not need to install any dependencies (they are included in the package).
Simply download the pre-compiled windows binary from here [windows binary](win_solution/iDemux_win10_64bit.zip).
Extract the zip file. To execute the tool press `windows+r`, enter `cmd`, `cd C:\\location_of_the_extracted_zip_file\bin` and execute `.\\idemuxCPP`


**From Source (distribution tar)**

To configure, compile and install execute the following commands on your command line:
```
    ./configure [--help for additional configuration options]
    make
    make install
```

**From Source (git)**
The installation from source requires additional tools and libraries:
* gnulib (autoconf, automake, etc.)

Generate the configure file with:
```
    autoreconf -i
```
Then proceed with `./configure` and `make` like in the previous [section](#from-source-(distribution-tar)).

**From Linux Package**
<table><thead><tr>
<th> Debian </th>
<th> Ubuntu </th>
<th> Fedora </th>
</tr></thead><tbody><tr>
<td style="vertical-align:top">
<details><summary>Debian_12</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Debian_12/i386/idemuxcpp_0.3.0-1_i386.deb"> idemuxcpp - 0.3.0 - 32 bit</a></p>
<p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Debian_12/amd64/idemuxcpp_0.3.0-1_amd64.deb"> idemuxcpp - 0.3.0 - 64 bit</a></p>
</details>
<details><summary>Debian_11</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Debian_11/i386/idemuxcpp_0.3.0-1_i386.deb"> idemuxcpp - 0.3.0 - 32 bit</a></p>
<p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Debian_11/amd64/idemuxcpp_0.3.0-1_amd64.deb"> idemuxcpp - 0.3.0 - 64 bit</a></p>
</details>
<details><summary>Debian_10</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Debian_10/i386/idemuxcpp_0.3.0-1_i386.deb"> idemuxcpp - 0.3.0 - 32 bit</a></p>
<p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Debian_10/amd64/idemuxcpp_0.3.0-1_amd64.deb"> idemuxcpp - 0.3.0 - 64 bit</a></p>
</details>
<details><summary>Debian_9.0</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Debian_9.0/i386/idemuxcpp_0.3.0-1_i386.deb"> idemuxcpp - 0.3.0 - 32 bit</a></p>
<p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Debian_9.0/amd64/idemuxcpp_0.3.0-1_amd64.deb"> idemuxcpp - 0.3.0 - 64 bit</a></p>
</details>
</td>
<td style="vertical-align:top">
<details><summary>xUbuntu_24.04</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/xUbuntu_24.04/amd64/idemuxcpp_0.3.0-1_amd64.deb"> idemuxcpp - 0.3.0 - 64 bit</a></p></details>
<details><summary>xUbuntu_23.10</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/xUbuntu_23.10/amd64/idemuxcpp_0.3.0-1_amd64.deb"> idemuxcpp - 0.3.0 - 64 bit</a></p></details>
<details><summary>xUbuntu_23.04</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/xUbuntu_23.04/amd64/idemuxcpp_0.3.0-1_amd64.deb"> idemuxcpp - 0.3.0 - 64 bit</a></p></details>
<details><summary>xUbuntu_22.10</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/xUbuntu_22.10/amd64/idemuxcpp_0.3.0-1_amd64.deb"> idemuxcpp - 0.3.0 - 64 bit</a></p></details>
<details><summary>xUbuntu_22.04</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/xUbuntu_22.04/amd64/idemuxcpp_0.3.0-1_amd64.deb"> idemuxcpp - 0.3.0 - 64 bit</a></p></details>
<details><summary>xUbuntu_20.04</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/xUbuntu_20.04/amd64/idemuxcpp_0.3.0-1_amd64.deb"> idemuxcpp - 0.3.0 - 64 bit</a></p></details>
</details><details><summary>xUbuntu_18.04</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/xUbuntu_18.04/amd64/idemuxcpp_0.3.0-1_amd64.deb"> idemuxcpp - 0.3.0 - 64 bit</a></p>
</details><details><summary>xUbuntu_16.04</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/xUbuntu_16.04/i386/idemuxcpp_0.3.0-1_i386.deb"> idemuxcpp - 0.3.0 - 32 bit</a></p>
<p><a href="https://download.opensuse.org/repositories/home:/Lexogen/xUbuntu_16.04/amd64/idemuxcpp_0.3.0-1_amd64.deb"> idemuxcpp - 0.3.0 - 64 bit</a></p></details>
</td>
<td style="vertical-align:top">
<details><summary>Fedora 40</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Fedora_40/x86_64/idemuxcpp-0.3.0-39.1.x86_64.rpm"> idemuxcpp - 0.3.0 - 64 bit</a></p></details>
<details><summary>Fedora 39</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Fedora_39/x86_64/idemuxcpp-0.3.0-39.1.x86_64.rpm"> idemuxcpp - 0.3.0 - 64 bit</a></p></details>
<details><summary>Fedora 38</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Fedora_38/x86_64/idemuxcpp-0.3.0-39.1.x86_64.rpm"> idemuxcpp - 0.3.0 - 64 bit</a></p></details>
<details><summary>Fedora 37</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Fedora_37/x86_64/idemuxcpp-0.3.0-39.1.x86_64.rpm"> idemuxcpp - 0.3.0 - 64 bit</a></p></details>
<details><summary>Fedora 36</summary><p><a href="https://download.opensuse.org/repositories/home:/Lexogen/Fedora_36/x86_64/idemuxcpp-0.3.0-39.1.x86_64.rpm"> idemuxcpp - 0.3.0 - 64 bit</a></p></details>
</td>
</tr></tbody></table>


on ubuntu you can install it for example with:
```
dpkg -i <idemuxcpp*.deb>
```

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
6. If your i5 has been sequenced as reverse complement, enter the reverse
   complement sequences in the sample sheet and use the ``--i5-rc`` option!
   If you are not sure, you can also use the ``--auto-detect`` option as alternative.


See [below](#sample-sheet-examples). for more showcases of sample/barcode combinations that are *allowed* or
*disallowed*.

## 3. Extract non-demultiplexed read data from a sequencing run
The read input files for idemux are non-demultiplexed read files which you can get by using demultiplexing software to extract reads from a sequencing run without demultiplexing by sample.  
You can use any demultiplexing software available to you, but the resulting read file(s) should contain all reads of the sequencing run you want to demultiplex with idemux.
Further, the reads should contain the read-out of the i7 + i5 barcode sequences in the read ID.

The following part of this section outlines how to use Illumina's bcl2fastq software to obtain the reads.

### Demultiplexing with bcl2fastq:
```
bcl2fastq -R /path/to/sequencing/run -o /path/to/output -l WARNING --no-lane-splitting --sample-sheet Illumina_EMPTY_SampleSheet.csv --barcode-mismatches 0 --mask-short-adapter-reads 10
```

This commands bcl2fastq to "demultiplex" the run at */path/to/sequencing/run* to the output directory */path/to/output*.
The content of the file *Illumina_EMPTY_SampleSheet.csv* has to match Illumina's format for the respective sequencer.

The following text is an example for the content of a SampleSheet for a Illumina Nextseq run:
```
[Header],,,,,,,
IEMFileVersion,4,,,,,,
Date,30.05.2017,,,,,,
Workflow,GenerateFASTQ,,,,,,
Application,NextSeq FASTQ Only,,,,,,
Assay,TruSeq RNA,,,,,,
Description,,,,,,,
Chemistry,Default,,,,,,
,,,,,,,
[Reads],,,,,,,
,,,,,,,
[Settings],,,,,,,
,,,,,,,
[Data],,,,,,,

Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
1,1,,,9999,AAAAAAAAAAAA,9999,AAAAAAAAAAAA,,
```

As you can see, no settings are specified and only one 'sample' was defined with a squence combination that is not likely to be close to any of the utilized barcode sequences.
**You have to adjust the length of the A\* stretches to the sequenced length of the i7/i5 barcodes!**
This specification is necessary to command bcl2fastq to write the i7+i5 sequence information in each read in the *Undetermined_S0_R1_001.fastq.gz* (*Undetermined_S0_R2_001.fastq.gz*) file(s)
The resulting reads in *Undetermined_S0_R1_001.fastq.gz* (*Undetermined_S0_R2_001.fastq.gz*) should follow this formatting style:

```
@NB502007:379:HM7H2BGXF:1:11101:19231:1159 1:N:0:TTAGGACGCAAA+GGGTCTGCCGAA
GCTCATCCATCTTTTTGAAAACTCTTCATACTCGTTAGATCGGAAGAG
+
AAAAAEEEAEEEEEEAEEEEEEEEEEEEEEEEEEEE/E/EEEEE/EEE
@NB502007:379:HM7H2BGXF:1:11101:17406:1159 1:N:0:AAGTAACAGCTT+AATCGTGGACGG
CACACCTCCGTTCACGACGCTCTTCCGATATAGATGTAACTGGAGGAA
+
AAAAAEEEEEAEE/EEEEEEEEEE/EEEEAEA/EEEEEEEEEEEEEEE
@NB502007:379:HM7H2BGXF:1:11101:18203:1159 1:N:0:CTGCCAACACGA+GCTGTGGTTCAT
GACATGTATACAGTCTACGGATGAACGTTTAGATCGGAAGAGCACACG
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEE
@NB502007:379:HM7H2BGXF:1:11101:7322:1159 1:N:0:TACATGGCCACT+ATGTTCCAGTGA
CTTGGTCACGCTACTGTACTCCAGCCAGGGCGACAGAGCAAGACCTAT
+
AAAAAEEEEEEEEEEEE/EEEEEEEEAEEEEEAEEEEEEEEEEEAEEE
...
 ```

## 4. Running idemuxCPP
Once you have installed the tool you can run it by typing ``idemuxCPP`` in the terminal.

idemuxCPP accepts the following arguments:
```
  -h, --help                    Print help and exit
  -V, --version                 Print version and exit

Required arguments:
  -1, --r1=STRING               Fastq.gz read file 1 (or .fastq file).
                                    (default='')
  -2, --r2=STRING               Fastq.gz read file 2 (required only in paired
                                  end mode).
                                    (default='')
  -o, --out=STRING              Where to write the output files.
                                    (default='./')
  -s, --sample-sheet=STRING     Input a csv file describing sample names and
                                  barcode combinations (i7, i5 and i1
                                  barcodes).
                                    (default='sample-sheet.csv')

Optional arguments:
  -b, --barcode-corrections=STRING
                                Outputs a csv file that contains the number of
                                  corrected barcodes
  -5, --i5-rc                   Should be set when the i5 barcode has been
                                  sequenced as reverse complement. Make sure to
                                  enter reverse complement sequences in the
                                  barcode file.  (default=off)
  -i, --i1-start=INT            Start position of the i1 index (1-based) on
                                  read 2.
                                    (default='11')
      --i1-read=INT             Read in which the i1 index should be corrected
                                  (1 or 2).
                                    (default='2')
  -q, --queue-size=INT          Queue size for reads that will be processed in
                                  one block.
                                    (default='4000000')
  -r, --reading-threads=INT     Number of threads used for reading gz files.
                                  Either 1 or 2 (one thread per input file is
                                  used).
                                    (default='2')
  -w, --writing-threads=INT     Number of threads used for writing gz files.
                                  Default is the number of processor cores.

  -p, --processing-threads=INT  Number of threads used for processing the error
                                  correction. Default is the number of
                                  processor cores.

  -d, --demux-only              Do a one on one mapping for the barcodes
                                  specified in the sample sheet. No error
                                  correction will be done. Barcodes that do not
                                  match are written to the undetermined reads
                                  file.  (default=off)
  -v, --verbose                 Verbose.
                                    (default=off)
```

Example commands:
```
    # demultiplexes read 1 and 2 into the folder 'demux'
    idemuxCPP --r1 read_1.fastq.gz --r2 read_2.fastq.gz --sample-sheet samples.csv --out demux

    # demultiplexing assuming the i1 barcode starts at the first base
    idemuxCPP --r1 read_1.fastq.gz --r2 read_2.fastq.gz --sample-sheet samples.csv --out demux --i1_start 1

    # demultiplexing assuming i5 is present as reverse complement in the fastq header
    # if he i5 has been sequenced as reverse complement use this option and provide
    # the reverse complement sequences in the sample sheet.
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

If you sequenced i5 as a reverse complement, make sure to fill in reverse complement
barcodes into the sample sheet and to use the ``--i5-rc`` parameter.

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

    # mixed indexing (if not ambiguous) (full i7 and sparse i5, i1)
    sample_name,i7,i5,i1
    sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
    sample_1,AAAATCCCAGTT,,AAAATCCCAGTT
    sample_2,GAAAATTTACGC,GCCCCTTTCAGA,GAAAATTTACGC
    sample_3,AAACTAACTGTC,,AAACTAACTGTC

    # mixed indexing (if not ambiguous) (no i7, sparse i5 & i1)
    sample_name,i7,i5,i1
    sample_0,,CCCCACTGAGTT,
    sample_1,,,AAAATCCCAGTT

    # mixed indexing (if not ambiguous) (sparse i7, full i5 & i1)
    sample_name,i7,i5,i1
    sample_0,,CCCCACTGAGTT,AAAACATGCGTT
    sample_1,AAAATCCCAGTT,CCCCTAAACGTT,AAAATCCCAGTT
    sample_2,,GCCCCTTTCAGA,GAAAATTTACGC
    sample_3,AAACTAACTGTC,CCCATCCATGTA,AAACTAACTGTC

    # additional parameter columns for i1_read and i1_start index (1-based).
    sample_name,i7,i5,i1,i1_read,i1_start
    sample_0,AANACATGCGTT,,TTTTAG,2,1
    sample_1,AANACATGCG,,AAAACATG,2,11
    sample_2,AANACA,,CACCCC,1,5
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

