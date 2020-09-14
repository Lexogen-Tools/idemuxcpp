#include <string>
#include <fstream>
#include <unordered_map>
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
	string csv_header = "sample\twritten_reads";
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

Writer::~Writer() {
}

