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

using namespace std;



int main(int argc, char **argv) {
	string read1_file = "";
	string read2_file = "";
	string outputdirectory = "";
	string sample_sheet_file = "";
	string barcode_corrections_file = "";
	bool i5_rc = false;
	int i1_start = 10; // zero based index
	size_t queue_size;
	int reading_threads;
	int writing_threads = -1;
	int processing_threads = -1;
	bool demux_only = false;
	//bool verbose = false;
	string relative_exepath = string(argv[0]);
	std::cout << relative_exepath << std::endl;

	struct idemuxCPP_args_info args_info;
	// check if there is any unparseable argument
	if (idemuxCPP_cmdline_parser(argc, argv, &args_info) != 0)
		exit(EXIT_FAILURE);

	if (args_info.r1_given) {
		read1_file = args_info.r1_arg;
	} else {
		fprintf(stderr, "Error: please enter a read1.fastq.gz file!\n");
		exit(EXIT_FAILURE);
	}
	if (args_info.r2_given) {
		read2_file = args_info.r2_arg;
	} else {
		fprintf(stderr, "Error: please enter a read2.fastq.gz file!\n");
		exit(EXIT_FAILURE);
	}
	if (args_info.out_given) {
		outputdirectory = args_info.out_arg;
	}
	if (args_info.sample_sheet_given) {
		sample_sheet_file = args_info.sample_sheet_arg;
	}

	queue_size = args_info.queue_size_arg;
	reading_threads = args_info.reading_threads_arg;

	if (args_info.writing_threads_given) {
		writing_threads = args_info.writing_threads_arg;
	}
	if (args_info.processing_threads_given) {
		processing_threads = args_info.processing_threads_arg;
	}
	if (args_info.barcode_corrections_given){
		barcode_corrections_file = string(args_info.barcode_corrections_arg);
	}
	demux_only = args_info.demux_only_flag;

	i5_rc = args_info.i5_rc_flag;
    //use a zero based starting position internally
	if (args_info.i1_start_arg < 1){
		fprintf(stderr, "Error: please enter a starting position >=1!\n");
		exit(EXIT_FAILURE);
	}
	i1_start = args_info.i1_start_arg -1;
	//verbose = args_info.verbose_flag;

	Parser p;
	vector<Barcode*> barcodes;
	unordered_map<string, string> *barcode_sample_map = p.parse_sample_sheet(
			sample_sheet_file, i5_rc, barcodes, relative_exepath, demux_only);

	// do things.
	demux_paired_end(barcode_sample_map, barcodes, read1_file, read2_file,
			i1_start, outputdirectory, p, queue_size, reading_threads, writing_threads, processing_threads, barcode_corrections_file);

	delete barcode_sample_map;
	for (size_t i = 0; i < barcodes.size(); i++)
		delete barcodes[i];
	idemuxCPP_cmdline_parser_free(&args_info);
	return EXIT_SUCCESS;
}

