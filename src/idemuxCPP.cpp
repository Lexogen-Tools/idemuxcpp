//============================================================================
// Name        : idemuxCPP.cpp
// Author      : GE
// Version     :
// Copyright   : (C) GE 2020
// Description : Demultiplex fastq.gz reads according to their indices.
//============================================================================

#include <zlib.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <cstring>
#include <algorithm>
#include <unordered_map>
#include <limits>
#include <cmath>

#include<idemuxCPP_cmdl.h>
#include "Parser.h"
#include "FileHandler.h"
#include "PairedReader.h"
#include "ZipFastqReader.h"
#include "ZipFastqWriter.h"

using namespace std;

/**
 * process mate pair

 Returns
 -------
 barcodes: ()
 */
string process_mate_pair(std::pair<fq_read*, fq_read*> &mate_pair,
		unordered_set<string> *i7_wanted, unordered_set<string> *i5_wanted,
		unordered_set<string> *i1_wanted,
		unordered_map<string, string> *map_i7,
		unordered_map<string, string> *map_i5,
		unordered_map<string, string> *map_i1, int i1_start, int i1_end,
		std::pair<fq_read*, fq_read*> &mate_pair_out) {

	string fastq_header(mate_pair.first->Seq_ID);
	std::pair<string, string> barcodes = Parser::parse_indices(fastq_header);
	// if not barcodes:
	string i7_bc, i5_bc;

	// when there are 2 barcodes in the fastq header the orientation is i7,i5
	string tmp_i7_bc = barcodes.first;
	string tmp_i5_bc = barcodes.second;
	if(i7_wanted->size() == 0)
		i7_bc = "";
	if(i5_wanted->size() == 0)
		i5_bc = "";

	i7_bc = map_i7->at(tmp_i7_bc);
	i5_bc = map_i5->at(tmp_i5_bc);

	string i1_bc = "";

	//if i7_bc in i7_wanted and i5_bc in i5_wanted:
	//if (it7 != i7_wanted->end() and it5 != i5_wanted->end()) {
	fq_read *r1c = new fq_read(*mate_pair.first);
	fq_read *r2c = new fq_read(*mate_pair.second);
	if(i1_wanted->size() > 0){
		i1_bc = mate_pair.second->Sequence.substr(i1_start, i1_end - i1_start);
		//i1_bc = mate_pair[1][1][i1_start:i1_end]
		auto iti1c = map_i1->find(i1_bc);
		if (iti1c != map_i1->end()) {
			string _i1_corrected = iti1c->second;
			//_i1_corrected = map_i1.get(i1_bc)

			//if _i1_corrected in i1_wanted:
			auto it1 = i1_wanted->find(_i1_corrected);
			if (it1 != i1_wanted->end()) {

				r1c->Seq_ID = string(mate_pair.first->Seq_ID) + "+" + i1_bc;
				//m1_hdr = f"{m1_hdr[:-1]}+{i1_bc}\n"
				r1c->Sequence = string(mate_pair.first->Sequence);
				r1c->Plus_ID = string(mate_pair.first->Plus_ID);
				r1c->QualityCode = string(mate_pair.first->QualityCode);

				r2c->Seq_ID = string(mate_pair.second->Seq_ID) + "+" + i1_bc;
				r2c->Sequence = string(mate_pair.second->Sequence).substr(0,
						i1_start)
						+ string(mate_pair.second->Sequence).substr(i1_end,
								mate_pair.second->Sequence.length()-i1_end);
				r2c->Plus_ID = string(mate_pair.second->Plus_ID);
				r2c->QualityCode = string(mate_pair.second->QualityCode).substr(
						0, i1_start)
						+ string(mate_pair.second->QualityCode).substr(i1_end,
								mate_pair.second->QualityCode.length() - i1_end);
				//m2_hdr = f"{m2_hdr[:-1]}+{i1_bc}\n"
				//m2_seq = f"{m2_seq[:i1_start]}{m2_seq[i1_end:]}"
				//m2_qcs = f"{m2_qcs[:i1_start]}{m2_qcs[i1_end:]}"
				i1_bc = _i1_corrected;
			}
		}
	}
	mate_pair_out.first = r1c;
	mate_pair_out.second = r2c;
	string barcodes_key = i7_bc + '\n' + i5_bc + '\n' + i1_bc;
	return barcodes_key;
}

void demux_paired_end(unordered_map<string, string> *barcode_sample_map,
		vector<Barcode*> &barcodes, string read1, string read2, int i1_start,
		string output_dir, Parser &parser) {
	// TODO: add documentation
	// TODO: add logging
	// load the maps that will be used for error correction. As the tool does not allow
	// different length we only need to load the used length
	Barcode *i7, *i5, *i1;
	i7 = barcodes[0];
	i5 = barcodes[1];
	i1 = barcodes[2];

	unordered_set<string> *i7_wanted = i7->used_codes();
	unordered_set<string> *i5_wanted = i5->used_codes();
	unordered_set<string> *i1_wanted = i1->used_codes();
	unordered_map<string, string> *map_i7 = i7->correction_map;
	unordered_map<string, string> *map_i5 = i5->correction_map;
	unordered_map<string, string> *map_i1 = i1->correction_map;
	int i1_end = i1_start + i1->length;

	// if None is in *_wanted no barcode has been specified
	bool has_i7 = !i7->empty();
	bool has_i5 = !i5->empty();
	bool has_i1 = !i1->empty();

	// before doing any processing check if the fastq file is okay.
	parser.peek_into_fastq_files(read1, read2, has_i7, has_i5, has_i1,
			i7->length, i5->length, i1_start, i1_end);

	//read_counter = Counter()
	std::cout << "Starting demultiplexing" << std::endl;
	// first we need to open the output files the reads should get sorted into
	FileHandler *file_handler = new FileHandler(*barcode_sample_map, output_dir,
			(size_t) pow(2, 30));

	//with FileHandler(barcode_sample_map, output_dir) as file_handler:
	// then we iterate over all the paired end reads
	PairedReader pr(read1, read2);
	std::vector<std::pair<fq_read*, fq_read*>> *pe_reads;
	bool reads_available = true;
	while (reads_available) {
		pe_reads = pr.next_reads(4000000);
		std::cout<< "n_reads " << pe_reads->size() << std::endl;
		if (pe_reads != NULL && pe_reads->size() > 0) {

			for (auto it = pe_reads->begin(); it != pe_reads->end(); it++) {
				std::pair<fq_read*, fq_read*> mate_pair = *it;
				std::pair<fq_read*, fq_read*> processed_mates;
				string s_barcodes = process_mate_pair(mate_pair, i7_wanted,
						i5_wanted, i1_wanted, map_i7, map_i5, map_i1,
						i1_start, i1_end, processed_mates);
				std::pair<ZipFastqWriter*, ZipFastqWriter*> *pair_writer =
						file_handler->get_file_handles(s_barcodes);
				//processed_mates
				if (pair_writer != NULL) {
					pair_writer->first->write_read(processed_mates.first);
					pair_writer->second->write_read(processed_mates.second);
				}
				delete processed_mates.first;
				delete processed_mates.second;
			}
		} else {
			reads_available = false;
		}
		if (pe_reads != NULL) {
			for (auto it = pe_reads->begin(); it != pe_reads->end(); it++) {
				delete it->first;
				delete it->second;
			}
			delete pe_reads;
		}

	}
	delete file_handler;
	delete i7_wanted;
	delete i5_wanted;
	delete i1_wanted;
	delete map_i7;
	delete map_i5;
	delete map_i1;
}


int main(int argc, char **argv) {
	string read1_file = "";
	string read2_file = "";
	string outputdirectory = "";
	string sample_sheet_file = "";
	bool i5_rc = false;
	int i1_start = 10;
	bool verbose = false;

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
	i5_rc = args_info.i5_rc_flag;
	i1_start = args_info.i1_start_arg;
	verbose = args_info.verbose_flag;

	Parser p;
	vector<Barcode*> barcodes;
	unordered_map<string, string> *barcode_sample_map = p.parse_sample_sheet(
			sample_sheet_file, i5_rc, barcodes);

	// do things.
	demux_paired_end(barcode_sample_map, barcodes, read1_file, read2_file,
			i1_start, outputdirectory, p);

	delete barcode_sample_map;
	for(int i = 0; i < barcodes.size(); i++) delete barcodes[i];
	idemuxCPP_cmdline_parser_free(&args_info);
	return EXIT_SUCCESS;
}

