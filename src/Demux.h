#ifndef SRC_DEMUX_H_
#define SRC_DEMUX_H_

#include "config.h"
#include <zlib.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <limits>
#include <cmath>
#include <chrono>

#include "Parser.h"
#include "FileHandler.h"
#include "FileHandlerSE.h"
#include "PairedReader.h"
#include "FastqReader.h"
#ifdef HAVE_LIBBAMTOOLS
#include "BamReader.h"
#endif
#include "BoostZipReader.h"
#include "ZipFastqWriter.h"
#include "Writer.h"
#include "Correction_Counter.h"
#ifdef _OPENMP //HAVE_OPENMP
#include "omp.h"
#endif

using namespace std;

typedef std::chrono::high_resolution_clock Clock;


string extract_i1(fq_read* r, size_t i1_start, size_t i1_end){
	return r->Sequence.substr(i1_start, i1_end - i1_start);
}
/**
 * remove the i1 interval from the sequence and the quality code and append it to the index.
 */
void cut_out_i1(fq_read* r, size_t i1_start, size_t i1_end){
	string i1_bc = r->Sequence.substr(i1_start, i1_end - i1_start);
	r->Seq_ID = string(r->Seq_ID).append("+").append(i1_bc);
	r->Sequence = string(r->Sequence).substr(0,i1_start).append(r->Sequence.substr(i1_end));
	//r->Plus_ID = string(r->Plus_ID);
	r->QualityCode = string(r->QualityCode).substr(0, i1_start).append(r->QualityCode.substr(i1_end));
}

/**
 * process mate pair

 Returns
 -------
 barcodes: ()
 */
string process_mate_pair(std::pair<fq_read*, fq_read*> &mate_pair,
		unordered_set<string> *i7_wanted, unordered_set<string> *i5_wanted,
		unordered_set<string> *i1_wanted, unordered_map<string, string> *map_i7,
		unordered_map<string, string> *map_i5,
		unordered_map<string, string> *map_i1, unordered_map<string, i1_info> &i7_i5_i1_info_map,
		std::pair<fq_read*, fq_read*> &mate_pair_out,
		Correction_Counter *counted_corrections_per_index,
                size_t max_length_i7 = SIZE_MAX,
                size_t max_length_i5 = SIZE_MAX) {
	string fastq_header(mate_pair.first->Seq_ID);
	std::pair<string, string> barcodes = Parser::parse_indices(fastq_header, max_length_i7, max_length_i5);
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

	bool corrected_i7, corrected_i5, corrected_i1;
	corrected_i7 = corrected_i5 = corrected_i1 = false;

	auto i7_bc_it = map_i7->find(tmp_i7_bc);
	if(i7_bc_it != map_i7->end()){
		i7_bc = i7_bc_it->second;
		if (counted_corrections_per_index != NULL && i7_bc.compare(tmp_i7_bc) != 0)
			corrected_i7 = true;
	}
	auto i5_bc_it = map_i5->find(tmp_i5_bc);
	if(i5_bc_it != map_i5->end()){
		i5_bc = i5_bc_it->second;
		if (counted_corrections_per_index != NULL && i5_bc.compare(tmp_i5_bc) != 0)
			corrected_i5 = true;
	}
	string i1_bc = "";


	//if i7_bc in i7_wanted and i5_bc in i5_wanted:
	//if (it7 != i7_wanted->end() and it5 != i5_wanted->end()) {
	fq_read *r1c = new fq_read(*mate_pair.first);
	fq_read *r2c = new fq_read(*mate_pair.second);
	/* correct i7 and i5 indices (write the original indices,
	   that could have a different length after parsing because of
	   the --restrict-barcode-length option)
	*/
	string new_seq_id = Parser::replace_indices(r1c->Seq_ID, tmp_i7_bc, tmp_i5_bc);
	r1c->Seq_ID = new_seq_id;
	new_seq_id = Parser::replace_indices(r2c->Seq_ID, tmp_i7_bc, tmp_i5_bc);
	r2c->Seq_ID = new_seq_id;

	if (i1_wanted->size() > 0) {
		string i7_i5_bc = i7_bc + "\n" + i5_bc;
		auto it_i1_info = i7_i5_i1_info_map.find(i7_i5_bc);
		if (it_i1_info != i7_i5_i1_info_map.end())
		{
			size_t i1_start = it_i1_info->second.start_index;
			size_t i1_end = it_i1_info->second.end_index;

			if(it_i1_info->second.read_index == 1)
				i1_bc = extract_i1(r1c, i1_start, i1_end);
			else
				i1_bc = extract_i1(r2c, i1_start, i1_end);

			auto iti1c = map_i1->find(i1_bc);
			if (iti1c != map_i1->end()) {
				string _i1_corrected = iti1c->second;
				//if _i1_corrected in i1_wanted:
				auto it1 = i1_wanted->find(_i1_corrected);
				if (it1 != i1_wanted->end()) {
					if(it_i1_info->second.read_index == 1){
						cut_out_i1(r1c, i1_start, i1_end);
						r2c->Seq_ID.append("+").append(i1_bc);
					}
					else{
						cut_out_i1(r2c, i1_start, i1_end);
						r1c->Seq_ID.append("+").append(i1_bc);
					}

					if (counted_corrections_per_index != NULL && i1_bc.compare(_i1_corrected) != 0)
						corrected_i1 = true;
					i1_bc = _i1_corrected;
				}
			}
		}
	}

	mate_pair_out.first = r1c;
	mate_pair_out.second = r2c;
	string barcodes_key = string("").append(i7_bc).append("\n").append(i5_bc).append("\n").append(
			i1_bc);

#pragma omp critical
	{
		if (counted_corrections_per_index != NULL)
			counted_corrections_per_index->count_correction(i7_bc,i5_bc,i1_bc,corrected_i7, corrected_i5, corrected_i1);
	}

	return barcodes_key;
}

string process_read(fq_read* read,
		unordered_set<string> *i7_wanted, unordered_set<string> *i5_wanted,
		unordered_set<string> *i1_wanted, unordered_map<string, string> *map_i7,
		unordered_map<string, string> *map_i5,
		unordered_map<string, string> *map_i1, unordered_map<string, i1_info> &i7_i5_i1_info_map,
		fq_read &read_out,
		Correction_Counter *counted_corrections_per_index,
                size_t max_length_i7,
                size_t max_length_i5) {
	string fastq_header(read->Seq_ID);
	std::pair<string, string> barcodes = Parser::parse_indices(fastq_header, max_length_i7, max_length_i5);
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

	bool corrected_i7, corrected_i5, corrected_i1;
	corrected_i7 = corrected_i5 = corrected_i1 = false;

	auto i7_bc_it = map_i7->find(tmp_i7_bc);
	if(i7_bc_it != map_i7->end()){
		i7_bc = i7_bc_it->second;
		if (counted_corrections_per_index != NULL && i7_bc.compare(tmp_i7_bc) != 0)
			corrected_i7 = true;
	}
	auto i5_bc_it = map_i5->find(tmp_i5_bc);
	if(i5_bc_it != map_i5->end()){
		i5_bc = i5_bc_it->second;
		if (counted_corrections_per_index != NULL && i5_bc.compare(tmp_i5_bc) != 0)
			corrected_i5 = true;
	}
	string i1_bc = "";

	read_out = fq_read(*read);

	/* correct i7 and i5 indices (write the original indices,
	   that could have a different length after parsing because of
	   the --restrict-barcode-length option)
	*/
	string new_seq_id = Parser::replace_indices(read->Seq_ID, tmp_i7_bc, tmp_i5_bc);
	read_out.Seq_ID = new_seq_id;

	if (i1_wanted->size() > 0) {
		string i7_i5_bc = i7_bc + "\n" + i5_bc;
		auto it_i1_info = i7_i5_i1_info_map.find(i7_i5_bc);
		if (it_i1_info != i7_i5_i1_info_map.end())
		{
			size_t i1_start = it_i1_info->second.start_index;
			size_t i1_end = it_i1_info->second.end_index;

			if(it_i1_info->second.read_index == 1)
				i1_bc = extract_i1(&read_out, i1_start, i1_end);

			auto iti1c = map_i1->find(i1_bc);
			if (iti1c != map_i1->end()) {
				string _i1_corrected = iti1c->second;
				//if _i1_corrected in i1_wanted:
				auto it1 = i1_wanted->find(_i1_corrected);
				if (it1 != i1_wanted->end()) {
					if(it_i1_info->second.read_index == 1){
						cut_out_i1(&read_out, i1_start, i1_end);
					}
					if (counted_corrections_per_index != NULL && i1_bc.compare(_i1_corrected) != 0)
						corrected_i1 = true;
					i1_bc = _i1_corrected;
				}
			}
		}
	}

	string barcodes_key = string("").append(i7_bc).append("\n").append(i5_bc).append("\n").append(
			i1_bc);
#pragma omp critical
	{
		if (counted_corrections_per_index != NULL)
			counted_corrections_per_index->count_correction(i7_bc,i5_bc,i1_bc,corrected_i7, corrected_i5, corrected_i1);
	}
	return barcodes_key;
}

size_t get_max_length(unordered_set<string> &barcodes){
    size_t max = 0;
    size_t bc_length = 0;
    for(auto it = barcodes.begin(); it != barcodes.end(); it++){
        bc_length = (*it).length();
        if(bc_length > max) { max = bc_length; }
    }
    return max;
}

void string_to_lower(string &out_string){
	transform(out_string.begin(), out_string.end(), out_string.begin(), ::tolower);
}

IFastqReader *init_reader_single_end(string reads_file){
	IFastqReader *reader = NULL;
	int last_index = reads_file.find_last_of('.');
	if(last_index > 0){
		string suffix = reads_file.substr(last_index);
		string suffix_lower = suffix;
		string_to_lower(suffix_lower);
		if(suffix_lower.compare(".bam") == 0){
#ifdef HAVE_LIBBAMTOOLS
			reader = new BamReader(reads_file);
#else
			throw runtime_error("Error: bam files are not supported with this compilation!");
#endif
		}
		else{
			if (suffix_lower.compare(".gz") == 0){
				reader = new BoostZipReader(reads_file);
			}
			else{
				// use fastq reader
				reader = new FastqReader(reads_file);
			}
		}
	}
	else{
		throw runtime_error("Error: no valid input file ending!");
	}
	return reader;
}

PairedReader *init_reader_paired_end(string read1_file, string read2_file){
	PairedReader *paired_reader = NULL;
	int last_index = read1_file.find_last_of('.');
	if(last_index > 0){
		string suffix = read1_file.substr(last_index);
		string suffix_lower = suffix;
		string_to_lower(suffix_lower);
		if(suffix_lower.compare(".bam") == 0){
			paired_reader = new PairedReader(read1_file);
		}
		else{
			paired_reader = new PairedReader(read1_file, read2_file);
		}
	}
	else{
		throw runtime_error("Error: no valid input file ending!");
	}
	return paired_reader;
}

void demux_paired_end(unordered_map<string, string> *barcode_sample_map,
		vector<Barcode*> &barcodes, string read1_file, string read2_file, unordered_map<string, i1_info> &i7_i5_i1_info_map,
		string output_dir, Parser &parser, size_t queue_size, int reading_threads, int writing_threads, int processing_threads, string barcode_corrections_file, 
                bool skip_check, bool restrict_barcode_length) {
	// load the maps that will be used for error correction. As the tool does not allow
	// different length we only need to load the used length
	Barcode *i7, *i5, *i1;
	i7 = barcodes[0];
	i5 = barcodes[1];
	i1 = barcodes[2];

	unordered_set<string> *i7_wanted = i7->used_codes();
	unordered_set<string> *i5_wanted = i5->used_codes();
	unordered_set<string> *i1_wanted = i1->used_codes();
	unordered_map<string, string> *map_i7 = &i7->Correction_map;
	unordered_map<string, string> *map_i5 = &i5->Correction_map;
	unordered_map<string, string> *map_i1 = &i1->Correction_map;

	// if None is in *_wanted no barcode has been specified
	bool has_i7 = !i7->empty();
	bool has_i5 = !i5->empty();
	bool has_i1 = !i1->empty();

	size_t max_length_i7 = SIZE_MAX;
	size_t max_length_i5 = SIZE_MAX;
	if (restrict_barcode_length){
		max_length_i7 = get_max_length(*i7_wanted);
		max_length_i5 = get_max_length(*i5_wanted);
	}
	// before doing any processing check if the fastq file is okay.
	if (!skip_check){
		PairedReader * pr = init_reader_paired_end(read1_file, read2_file);
		parser.peek_into_fastq_files(*pr, has_i7, has_i5, has_i1,
				i7->Lengths, i5->Lengths, i7_i5_i1_info_map, max_length_i7, max_length_i5);
		delete pr;
	}
	Correction_Counter* counted_corrections_per_index = NULL;
	if(barcode_corrections_file.compare("") != 0){
		counted_corrections_per_index = new Correction_Counter(i7_wanted,i5_wanted,i1_wanted);
	}

	//read_counter = Counter()
	std::cout << "Starting demultiplexing" << std::endl;
	PairedReader * pr = init_reader_paired_end(read1_file, read2_file);
	// first we need to open the output files the reads should get sorted into
	FileHandler *file_handler = new FileHandler(*barcode_sample_map, output_dir,
			(size_t) pow(2, 30));

	// then we iterate over all the paired end reads
	std::vector<std::pair<fq_read*, fq_read*>> *pe_reads;
	bool reads_available = true;
	std::unordered_map<std::string, size_t> read_counter;
#ifdef _OPENMP
	int nproc_writing = omp_get_num_procs();
	if (writing_threads > 0)
		nproc_writing = writing_threads;
	int nproc_processing = omp_get_num_procs();
	if (processing_threads > 0)
		nproc_processing = processing_threads;
#endif
	while (reads_available) {
		auto t1 = Clock::now();
		pe_reads = pr->next_reads2(queue_size, reading_threads);
		auto t2 = Clock::now();

		std::cout << "n_reads " << pe_reads->size() << std::endl;
		std::cout << "time reading (" << pe_reads->size() << " reads): " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds" << std::endl;

		t1 = Clock::now();
		unordered_map<std::pair<ZipFastqWriter*, ZipFastqWriter*>*,
				vector<std::pair<fq_read*, fq_read*>>> map_pairs;
		if (pe_reads != NULL && pe_reads->size() > 0) {
			int max_reads = (int)pe_reads->size();
#pragma omp parallel for num_threads(nproc_processing) shared(pe_reads, i7_wanted, i5_wanted, i1_wanted,map_i7, map_i5, map_i1, i7_i5_i1_info_map, file_handler, map_pairs, counted_corrections_per_index)
			for (int i = 0; i < max_reads; i++){
				std::pair<fq_read*, fq_read*> mate_pair = pe_reads->at(i);
				std::pair<fq_read*, fq_read*> processed_mates;
				string s_barcodes = process_mate_pair(mate_pair, i7_wanted, i5_wanted,
						i1_wanted, map_i7, map_i5, map_i1, i7_i5_i1_info_map,
						processed_mates, counted_corrections_per_index, max_length_i7, max_length_i5);
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
					s1 = itread->first->to_string();
					s2 = itread->second->to_string();
				        key->first->write(s1.c_str(), s1.length());
				        key->second->write(s2.c_str(), s2.length());
					delete itread->first;
					delete itread->second;
				}
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

	if(barcode_corrections_file.compare("") != 0){
		counted_corrections_per_index->write_correction(barcode_corrections_file, *barcode_sample_map);
		delete counted_corrections_per_index;
	}
	Writer w;
	w.write_summary(read_counter, output_dir);

	delete pr;
	delete file_handler;
	delete i7_wanted;
	delete i5_wanted;
	delete i1_wanted;
}

void demux_single_end(unordered_map<string, string> *barcode_sample_map,
		vector<Barcode*> &barcodes, string reads_file, unordered_map<string, i1_info> &i7_i5_i1_info_map,
		string output_dir, Parser &parser, size_t queue_size, int reading_threads, int writing_threads, int processing_threads, string barcode_corrections_file, 
                bool skip_check, bool restrict_barcode_length) {
	// load the maps that will be used for error correction. As the tool does not allow
	// different length we only need to load the used length
	Barcode *i7, *i5, *i1;
	i7 = barcodes[0];
	i5 = barcodes[1];
	i1 = barcodes[2];

	unordered_set<string> *i7_wanted = i7->used_codes();
	unordered_set<string> *i5_wanted = i5->used_codes();
	unordered_set<string> *i1_wanted = i1->used_codes();
	unordered_map<string, string> *map_i7 = &i7->Correction_map;
	unordered_map<string, string> *map_i5 = &i5->Correction_map;
	unordered_map<string, string> *map_i1 = &i1->Correction_map;

	// if None is in *_wanted no barcode has been specified
	bool has_i7 = !i7->empty();
	bool has_i5 = !i5->empty();
	bool has_i1 = !i1->empty();

        size_t max_length_i7 = SIZE_MAX;
        size_t max_length_i5 = SIZE_MAX;
        if (restrict_barcode_length){
            max_length_i7 = get_max_length(*i7_wanted);
            max_length_i5 = get_max_length(*i5_wanted);
        }

	// before doing any processing check if the fastq file is okay.
        if (!skip_check){
        	IFastqReader *reader = init_reader_single_end(reads_file);
        	parser.peek_into_fastq_file(reader, has_i7, has_i5, has_i1,
			    i7->Lengths, i5->Lengths, i7_i5_i1_info_map, max_length_i7, max_length_i5);
        	delete reader;
        }
	Correction_Counter* counted_corrections_per_index = NULL;
	if(barcode_corrections_file.compare("") != 0){
		counted_corrections_per_index = new Correction_Counter(i7_wanted,i5_wanted,i1_wanted);
	}

	//read_counter = Counter()
	std::cout << "Starting demultiplexing" << std::endl;
	IFastqReader *reader = init_reader_single_end(reads_file);
	// first we need to open the output files the reads should get sorted into
	FileHandlerSE *file_handler = new FileHandlerSE(*barcode_sample_map, output_dir,
			(size_t) pow(2, 30));

	// then we iterate over all the paired end reads
	bool reads_available = true;
	std::unordered_map<std::string, size_t> read_counter;
#ifdef _OPENMP
	int nproc_writing = omp_get_num_procs();
	if (writing_threads > 0)
		nproc_writing = writing_threads;
	int nproc_processing = omp_get_num_procs();
	if (processing_threads > 0)
		nproc_processing = processing_threads;
#endif
	while (reads_available) {
		auto t1 = Clock::now();
		std::vector<fq_read*> se_reads;
		for(size_t i = 0; i < queue_size; i++){
			fq_read* r = reader->next_read();
			if(r){
				se_reads.push_back(r);
			}
			else{
				break;
			}
		}
		auto t2 = Clock::now();
		std::cout << "n_reads " << se_reads.size() << std::endl;
		std::cout << "time reading (" << se_reads.size() << " reads): " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds" << std::endl;

		t1 = Clock::now();
		unordered_map<ZipFastqWriter*,
				vector<fq_read*>> map_reads;
		if (se_reads.size() > 0) {
			int max_reads = (int)se_reads.size();
#pragma omp parallel for num_threads(nproc_processing) shared(se_reads, i7_wanted, i5_wanted, i1_wanted,map_i7, map_i5, map_i1, i7_i5_i1_info_map, file_handler, map_reads, counted_corrections_per_index)
			for (int i = 0; i < max_reads; i++){
				fq_read* r = se_reads[i];
				fq_read* processed_read = new fq_read();
				string s_barcodes = process_read(r, i7_wanted, i5_wanted,
						i1_wanted, map_i7, map_i5, map_i1, i7_i5_i1_info_map,
						*processed_read, counted_corrections_per_index, max_length_i7, max_length_i5);
				ZipFastqWriter *read_writer = file_handler->get_file_handles(s_barcodes);
#pragma omp critical
				map_reads[read_writer].push_back(processed_read);
			}

			vector<ZipFastqWriter*> keys;
			for (auto it = map_reads.begin(); it != map_reads.end(); it++) {
				keys.push_back(it->first);
			}
			t2 = Clock::now();
			std::cout << "time mapping: " << std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count() << " seconds" << std::endl;

			t1 = Clock::now();
#pragma omp parallel for num_threads(nproc_writing) shared(read_counter, keys, map_reads, file_handler)
			for (int i = 0; i < (int)keys.size(); i++) {
				ZipFastqWriter *key = keys[i];
				vector<fq_read*> val = map_reads[key];
				string s1 = "";
				string s2 = "";
				for (auto itread = val.begin(); itread != val.end(); itread++) {
					fq_read* r = (*itread);
					s1 += r->to_string();
					delete r;
				}
				key->write(s1.c_str(), s1.length());
			}
			// count the number of reads in an separate loop (maybe faster without locks)
			for (int i = 0; i < (int)keys.size(); i++) {
				ZipFastqWriter *key = keys[i];
				string sample_name = file_handler->get_sample_name(key);
				vector<fq_read*> val = map_reads[key];
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

		for (auto it = se_reads.begin(); it != se_reads.end(); it++) {
			delete (*it);
		}
	}

	if(barcode_corrections_file.compare("") != 0){
               counted_corrections_per_index->write_correction(barcode_corrections_file, *barcode_sample_map);
		delete counted_corrections_per_index;
	}
	Writer w;
	w.write_summary(read_counter, output_dir);

	delete reader;
	delete file_handler;
	delete i7_wanted;
	delete i5_wanted;
	delete i1_wanted;
}


#endif /* SRC_DEMUX_H_ */
