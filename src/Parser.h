/*
 * Parser.h
 *
 *  Created on: Aug 15, 2020
 *      Author: quin
 */

#ifndef PARSER_H_
#define PARSER_H_

#include <string>
#include <vector>
#include "ZipFastqReader.h"
#include "Barcode.h"

using namespace std;

class Parser {
public:
	Parser();
	unordered_map<string, string>* parse_sample_sheet(string sample_sheet,
			bool i5_rc, vector<Barcode*> &barcodes_out);
	string reverse_complement(string sequence);
	bool has_valid_barcode_combinations(Barcode &i7, Barcode &i5, Barcode &i1);

	void get_pe_fastq(string fq_gz_1, string fq_gz_2);
	void fastq_lines_to_reads(string fastq_lines);
	unordered_map<string, string>* get_map_from_resource(string package,
			string resource);
	Barcode* load_correction_map(Barcode &barcode);
	void peek_into_fastq_files(string fq_gz_1, string fq_gz_2, bool has_i7,
			bool has_i5, bool has_i1, int i7_length, int i5_length,
			int i1_start, int i1_end);
	void check_mate_pair(std::pair<fq_read*, fq_read*> mate_pair, bool has_i7,
			bool has_i5, bool has_i1, int i7_length, int i5_length,
			int i1_start, int i1_end);
	void check_mate2_length(fq_read *mate2, int i1_start, int i1_end);
	void check_fastq_headers(std::pair<fq_read*, fq_read*> mate_pair,
			bool has_i7, bool has_i5, bool i7_length, int i5_length);

	static std::pair<string, string> parse_indices(string input) {
		size_t index_colon = input.find(":");
		size_t index_plus1 = input.find("+");
		static std::pair<string, string> bcs_mate1 = std::pair<string, string>(
				input.substr(index_colon + 1, index_plus1),
				input.substr(index_plus1 + 1, input.length()));
		return bcs_mate1;
	}

	virtual ~Parser();
};

#endif /* PARSER_H_ */

