#ifndef SRC_DEMUX_H_
#define SRC_DEMUX_H_

#include "config.h"
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
#include <chrono>

#include "Parser.h"
#include "FileHandler.h"
#include "PairedReader.h"
#include "ZipFastqReader.h"
#include "ZipFastqWriter.h"
#include "Writer.h"
#ifdef HAVE_OPENMP
#include "omp.h"
#endif

using namespace std;

typedef std::chrono::high_resolution_clock Clock;

/**
 * process mate pair

 Returns
 -------
 barcodes: ()
 */
inline
string process_mate_pair(std::pair<fq_read*, fq_read*> &mate_pair,
		unordered_set<string> *i7_wanted, unordered_set<string> *i5_wanted,
		unordered_set<string> *i1_wanted, unordered_map<string, string> *map_i7,
		unordered_map<string, string> *map_i5,
		unordered_map<string, string> *map_i1, int i1_start, int i1_end,
		std::pair<fq_read*, fq_read*> &mate_pair_out) {

	string fastq_header(mate_pair.first->Seq_ID);
	std::pair<string, string> barcodes = Parser::parse_indices(fastq_header);
	// if not barcodes:
	string i7_bc = "";
	string i5_bc = "";

	// when there are 2 barcodes in the fastq header the orientation is i7,i5
	string tmp_i7_bc = barcodes.first;
	string tmp_i5_bc = barcodes.second;
	if (i7_wanted->size() == 0)
		i7_bc = "";
	if (i5_wanted->size() == 0)
		i5_bc = "";

	auto i7_bc_it = map_i7->find(tmp_i7_bc);
	if(i7_bc_it != map_i7->end())
		i7_bc = i7_bc_it->second;
	auto i5_bc_it = map_i5->find(tmp_i5_bc);
	if(i5_bc_it != map_i5->end())
		i5_bc = i5_bc_it->second;
	string i1_bc = "";

	//if i7_bc in i7_wanted and i5_bc in i5_wanted:
	//if (it7 != i7_wanted->end() and it5 != i5_wanted->end()) {
	fq_read *r1c = new fq_read(*mate_pair.first);
	fq_read *r2c = new fq_read(*mate_pair.second);
	if (i1_wanted->size() > 0) {
		i1_bc = mate_pair.second->Sequence.substr(i1_start, i1_end - i1_start);
		//i1_bc = mate_pair[1][1][i1_start:i1_end]
		auto iti1c = map_i1->find(i1_bc);
		if (iti1c != map_i1->end()) {
			string _i1_corrected = iti1c->second;
			//_i1_corrected = map_i1.get(i1_bc)

			//if _i1_corrected in i1_wanted:
			auto it1 = i1_wanted->find(_i1_corrected);
			if (it1 != i1_wanted->end()) {

				r1c->Seq_ID =
						string(mate_pair.first->Seq_ID).append("+").append(
								i1_bc);
				//m1_hdr = f"{m1_hdr[:-1]}+{i1_bc}\n"
				r1c->Sequence = string(mate_pair.first->Sequence);
				r1c->Plus_ID = string(mate_pair.first->Plus_ID);
				r1c->QualityCode = string(mate_pair.first->QualityCode);

				r2c->Seq_ID =
						string(mate_pair.second->Seq_ID).append("+").append(
								i1_bc);
				r2c->Sequence = string(mate_pair.second->Sequence).substr(0,
						i1_start).append(
						mate_pair.second->Sequence.substr(i1_end,
								mate_pair.second->Sequence.length() - i1_end));
				r2c->Plus_ID = string(mate_pair.second->Plus_ID);
				r2c->QualityCode = string(mate_pair.second->QualityCode).substr(
						0, i1_start).append(
						mate_pair.second->QualityCode.substr(i1_end,
								mate_pair.second->QualityCode.length()
										- i1_end));
				//m2_hdr = f"{m2_hdr[:-1]}+{i1_bc}\n"
				//m2_seq = f"{m2_seq[:i1_start]}{m2_seq[i1_end:]}"
				//m2_qcs = f"{m2_qcs[:i1_start]}{m2_qcs[i1_end:]}"
				i1_bc = _i1_corrected;
			}
		}
	}
	mate_pair_out.first = r1c;
	mate_pair_out.second = r2c;
	string barcodes_key = i7_bc.append("\n").append(i5_bc).append("\n").append(
			i1_bc);
	return barcodes_key;
}

inline
void demux_paired_end(unordered_map<string, string> *barcode_sample_map,
		vector<Barcode*> &barcodes, string read1, string read2, int i1_start,
		string output_dir, Parser &parser, size_t queue_size, int reading_threads, int writing_threads, int processing_threads) {
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

	// then we iterate over all the paired end reads
	PairedReader pr(read1, read2);
	std::vector<std::pair<fq_read*, fq_read*>> *pe_reads;
	bool reads_available = true;
	std::unordered_map<std::string, size_t> read_counter;
#ifdef HAVE_OPENMP
	int nproc_writing = omp_get_num_procs();
	if (writing_threads > 0)
		nproc_writing = writing_threads;
	int nproc_processing = omp_get_num_procs();
	if (processing_threads > 0)
		nproc_processing = processing_threads;
#endif
	while (reads_available) {
		auto t1 = Clock::now();
		pe_reads = pr.next_reads2(queue_size, reading_threads);
		auto t2 = Clock::now();

		std::cout << "n_reads " << pe_reads->size() << std::endl;
		std::cout << "time reading (" << pe_reads->size() << " reads): " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds" << std::endl;

		t1 = Clock::now();
		unordered_map<std::pair<ZipFastqWriter*, ZipFastqWriter*>*,
				vector<std::pair<fq_read*, fq_read*>>> map_pairs;
		if (pe_reads != NULL && pe_reads->size() > 0) {
#pragma omp parallel for num_threads(nproc_processing) shared(pe_reads, i7_wanted, i5_wanted, i1_wanted,map_i7, map_i5, map_i1, file_handler, map_pairs)
			for (int i = 0; i < (int)pe_reads->size(); i++){
				std::pair<fq_read*, fq_read*> mate_pair = pe_reads->at(i);
				std::pair<fq_read*, fq_read*> processed_mates;
				string s_barcodes = process_mate_pair(mate_pair, i7_wanted, i5_wanted,
						i1_wanted, map_i7, map_i5, map_i1, i1_start, i1_end,
						processed_mates);
				std::pair<ZipFastqWriter*, ZipFastqWriter*> *pair_writer = file_handler->get_file_handles(s_barcodes);
#pragma omp critical
				map_pairs[pair_writer].push_back(processed_mates);
			}

			vector<std::pair<ZipFastqWriter*, ZipFastqWriter*>*> keys;
			for (auto it = map_pairs.begin(); it != map_pairs.end(); it++) {
				keys.push_back(it->first);
			}
			t2 = Clock::now();
			std::cout << "time mapping: " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds" << std::endl;

			t1 = Clock::now();
#pragma omp parallel for num_threads(nproc_writing) shared(read_counter, keys, map_pairs, file_handler)
			for (int i = 0; i < (int)keys.size(); i++) {
				std::pair<ZipFastqWriter*, ZipFastqWriter*> *key = keys[i];
				vector<std::pair<fq_read*, fq_read*>> val = map_pairs[key];
				string s1 = "";
				string s2 = "";
				for (auto itread = val.begin(); itread != val.end(); itread++) {
					s1 += itread->first->to_string();
					s2 += itread->second->to_string();
					delete itread->first;
					delete itread->second;
				}
				key->first->write(s1.c_str(), s1.length());
				key->second->write(s2.c_str(), s2.length());
			}
			// count the number of reads in an separate loop (maybe faster without locks)
			for (int i = 0; i < (int)keys.size(); i++) {
				std::pair<ZipFastqWriter*, ZipFastqWriter*> *key = keys[i];
				string sample_name = file_handler->get_sample_name(key);
				vector<std::pair<fq_read*, fq_read*>> val = map_pairs[key];
				auto it_rc = read_counter.find(sample_name);
				if (it_rc == read_counter.end())
					read_counter[sample_name] = val.size();
				else
					it_rc->second += val.size();
			}
			t2 = Clock::now();
			std::cout << "time writing: " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds" << std::endl;

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

	Writer w;
	w.write_summary(read_counter, output_dir);

	delete file_handler;
	delete i7_wanted;
	delete i5_wanted;
	delete i1_wanted;
}



#endif /* SRC_DEMUX_H_ */
