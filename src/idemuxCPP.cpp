//============================================================================
// Name        : idemuxCPP.cpp
// Author      : GE
// Version     :
// Copyright   : (C) GE 2020
// Description : Demultiplex fastq.gz reads according to their indices.
//============================================================================

#include <iostream>
#include <cstdlib>
#include <stdio.h>

#include "idemuxCPP_cmdl.h"
#include "Parser.h"
#include "Barcode.h"
#include "Demux.h"
#include "helper.h"
#include "config.h"
#include <thread>

using namespace std;


int
main(int  argc,
     char **argv)
{
  string                      read1_file                = "";
  string                      read2_file                = "";
  string                      outputdirectory           = "";
  string                      sample_sheet_file         = "";
  string                      barcode_corrections_file  = "";
  string                      correction_maps_path      = "";
  bool                        i5_rc                     = false;
  bool                        i7_rc                     = false;
  bool                        auto_detect               = false;
  bool                        single_end_mode           = false;
  int                         default_i1_read           = 2;  //read in which the i1 index should be corrected (1 or 2).
  int                         default_i1_start          = 10; // zero based index
  size_t                      queue_size;
  size_t                      reading_threads;
  size_t                      n                       = (size_t)std::thread::hardware_concurrency();
  size_t                      writing_threads         = n;
  size_t                      processing_threads      = n;
  bool                        demux_only              = false;
  bool                        skip_check              = false;
  bool                        restrict_barcode_length = false;
  double                      queue_buffer_gb;
  size_t                      gzip_block_size_bytes = 1024 * 1024;
  bool                        verbose               = false;

  string                      relative_exepath = string(argv[0]);

  std::cout << relative_exepath << std::endl;

  struct idemuxCPP_args_info  args_info;

  // check if there is any unparseable argument
  if (idemuxCPP_cmdline_parser(argc, argv, &args_info) != 0)
    exit(EXIT_FAILURE);

  skip_check = args_info.skip_check_flag;

  single_end_mode = !args_info.paired_flag;
  if (args_info.bam_given) {
    read1_file = args_info.bam_arg;
  } else {
    if (args_info.r1_given) {
      read1_file = args_info.r1_arg;
    } else {
      fprintf(stderr,
              "Info: no read1 fastq file given! Trying to read from stdin without format check.\n");
      skip_check = true;
    }

    if (args_info.r2_given) {
      read2_file      = args_info.r2_arg;
      single_end_mode = false;
    } else {
      single_end_mode = true;
      default_i1_read = 1;
      if (!args_info.paired_flag)
        std::cout <<
        "Info: no read2.fastq.gz file given! Single end demultiplexing will be used." << std::endl;
    }
  }

  if (args_info.out_given)
    outputdirectory = args_info.out_arg;

  if (args_info.sample_sheet_given)
    sample_sheet_file = args_info.sample_sheet_arg;

  queue_size      = args_info.queue_size_arg;
  queue_buffer_gb = args_info.queue_buffer_gb_arg;
  reading_threads = (size_t)args_info.reading_threads_arg;

  if (args_info.writing_threads_given)
    writing_threads = (size_t)args_info.writing_threads_arg;

  if (args_info.processing_threads_given)
    processing_threads = (size_t)args_info.processing_threads_arg;

  if (args_info.barcode_corrections_given)
    barcode_corrections_file = string(args_info.barcode_corrections_arg);

  if (args_info.correction_map_prefix_given)
    correction_maps_path = string(args_info.correction_map_prefix_arg);

  demux_only = args_info.demux_only_flag;

  restrict_barcode_length = args_info.restrict_barcode_length_flag;

  i5_rc       = args_info.i5_rc_flag;
  i7_rc       = args_info.i7_rc_flag;
  auto_detect = args_info.auto_detect_flag;


  if (args_info.i1_start_arg < 1) {
    fprintf(stderr, "Error: please enter a starting position >=1!\n");
    exit(EXIT_FAILURE);
  }

  //use a zero based starting position internally
  default_i1_start = args_info.i1_start_arg - 1;

  if (args_info.i1_read_arg != 1 && args_info.i1_read_arg != 2) {
    fprintf(stderr, "Error: the i1_read parameter must be 1 or 2!\n");
    exit(EXIT_FAILURE);
  }

  default_i1_read       = args_info.i1_read_arg;
  gzip_block_size_bytes = args_info.gzip_block_size_arg;
  verbose               = args_info.verbose_flag;

  Parser                          p;
  vector<Barcode *>               barcodes;
  unordered_map<string, i1_info>  i7_i5_i1_info_map;
  unordered_map<string, string>   *barcode_sample_map = p.parse_sample_sheet(
    sample_sheet_file,
    i5_rc,
    i7_rc,
    auto_detect,
    barcodes,
    i7_i5_i1_info_map,
    relative_exepath,
    correction_maps_path,
    demux_only,
    default_i1_read,
    default_i1_start,
    single_end_mode);

  // do things.
  demux(barcode_sample_map, barcodes, read1_file, read2_file,
        i7_i5_i1_info_map, outputdirectory, p, queue_size, reading_threads, writing_threads,
        processing_threads, barcode_corrections_file, skip_check, restrict_barcode_length,
        queue_buffer_gb, gzip_block_size_bytes, single_end_mode, verbose);

  delete barcode_sample_map;
  for (size_t i = 0; i < barcodes.size(); i++)
    delete barcodes[i];
  idemuxCPP_cmdline_parser_free(&args_info);
  return EXIT_SUCCESS;
}
