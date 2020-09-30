#include <string>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include "Writer.h"
#include "Parser.h"

using namespace std;

Writer::Writer() {
}

/**
 * Writes a counter dictionary as tsv file with the name error_correction_stats.tsv.
 *
 * @param counter(dict): A dictionary barcodes <sample name, #corrected reads>.
 * @param output_dir (string): The path the file should be written to.
 *
 */
void Writer::write_summary(std::unordered_map<string, size_t> &counter,
		string output_dir) {
	string output_file = string("demultipexing_stats.tsv");
	if (output_dir.length() > 0)
		output_file = output_dir + PATH_SEP + string("demultipexing_stats.tsv");
	ofstream csvfile(output_file, std::ofstream::out);
	string csv_header = "sample_name\twritten_reads";
	csvfile << csv_header << std::endl;
	string sample_filename;
	size_t counts;
	for (auto it = counter.begin(); it != counter.end(); it++) {
		sample_filename = it->first;
		counts = it->second;
		csvfile << sample_filename << "\t" << counts << std::endl;
	}
	fprintf(stdout, "Run complete! Summary statistics saved to %s\n",
			output_file.c_str());
}

struct {
	bool operator()(const pair<string, size_t> &a, const pair<string, size_t> &b) const{
		if (a.second == b.second)
			return a.first.compare(b.first);
		return a.second < b.second;
	}
} compare_barcode_count;

void Writer::write_barcode_summary(
		correction_counter *counted_corrections_per_index,
		std::string barcode_corrections_file) {
	ofstream csvfile(barcode_corrections_file, std::ofstream::out);

	unordered_map<string, size_t>* count_map;

	string csv_header = "barcode_type\tbarcode\tcount_corrected";
	csvfile << csv_header << std::endl;
	string barcode_type;
	string barcode;
	size_t counts;
	for (int i = 0; i < 3; i++){
		unordered_map<string, size_t>* count_map;
		switch(i){
		case 0:
			barcode_type = "i1";
			count_map = &counted_corrections_per_index->count_i1;
			break;
		case 1:
			barcode_type = "i5";
			count_map = &counted_corrections_per_index->count_i5;
			break;
		case 2:
			barcode_type = "i7";
			count_map = &counted_corrections_per_index->count_i7;
			break;
		default:
			break;
		}
		vector<pair<string, size_t>> sorted;
		for (unordered_map<string, size_t>::iterator it = count_map->begin(); it != count_map->end(); it++)
			sorted.push_back(pair<string, size_t>(it->first, it->second));
		sort(sorted.begin(),sorted.end(),compare_barcode_count);
		for (auto it = sorted.begin(); it != sorted.end(); it++) {
			barcode = it->first;
			counts = it->second;
			csvfile << barcode_type << "\t"<< barcode << "\t" << counts << std::endl;
		}
	}
	fprintf(stdout, "Barcode corrections summary statistics saved to %s\n",
			barcode_corrections_file.c_str());
}


Writer::~Writer() {
}

