#ifndef PARSER_H_
#define PARSER_H_

#include <string>
#include <vector>
#include <fstream>
#include "FastqReader.h"
#include "Barcode.h"
#include "helper.h"

#ifdef _WIN32
 #include "direct.h"
 #define PATH_SEP '\\'
 #define GETCWD _getcwd
 #define CHDIR _chdir
#else
 #include "unistd.h"
 #define PATH_SEP '/'
 #define GETCWD getcwd
 #define CHDIR chdir
#endif
#include <sys/stat.h>

using namespace std;


class Parser {
public:
	Parser();
	std::vector<std::vector<std::string>> readCSV(std::istream &in);
	unordered_map<string, string>* parse_sample_sheet(string sample_sheet,
			bool i5_rc, vector<Barcode*> &barcodes_out, string relative_exepath, bool demux_only);
	string reverse_complement(string sequence);
	bool has_valid_barcode_combinations(Barcode &i7, Barcode &i5, Barcode &i1);
	void fastq_lines_to_reads(string fastq_lines);
	void peek_into_fastq_files(string fq_gz_1, string fq_gz_2, bool has_i7,
			bool has_i5, bool has_i1, int i7_length, int i5_length,
			int i1_start, int i1_end);
	void check_mate_pair(std::pair<fq_read*, fq_read*> mate_pair, bool has_i7,
			bool has_i5, bool has_i1, int i7_length, int i5_length,
			int i1_start, int i1_end);
	void check_mate2_length(fq_read *mate2, int i1_start, int i1_end);
	void check_fastq_headers(std::pair<fq_read*, fq_read*> mate_pair,
			bool has_i7, bool has_i5, int i7_length, int i5_length);

	static std::pair<string, string> parse_indices(string input) {
		size_t index_colon = input.find_last_of(':');
		size_t index_plus1 = input.find('+');
		string code_i7 = "";
		string code_i5 = "";
		if(index_colon!=std::string::npos)
			code_i7 = input.substr(index_colon + 1, index_plus1-index_colon-1);
		if(index_plus1!=std::string::npos)
			code_i5 = input.substr(index_plus1 + 1, input.length()-index_plus1-1);
		std::pair<string, string> bcs_mate1 = std::pair<string, string>(
				code_i7,
				code_i5);
		return bcs_mate1;
	}

	static
	bool PathExists(const std::string &s)
	{
	  struct stat buffer;
	  return (stat (s.c_str(), &buffer) == 0);
	}

	/*
	 * load error correction map from misc folder.
	 */
	static
	unordered_map<string, string>* get_map_from_resource(string package,
			string resource) {
		printf("Loading error correction map from %s\n", resource.c_str());
		unordered_map<string, string> *mapping =
				new unordered_map<string, string>();

		ifstream dataFile;
		string filepath = package + PATH_SEP + resource;
		dataFile.open(filepath);
		fprintf(stdout, "checking file for valid barcodes: %s\n", filepath.c_str());
		int tabidx;
		int return_idx; //for windows files
		string v1, v2;
		if (dataFile.fail()) {
			string message = string_format(
					"Error: failed to open the correction map file! Please copy it "
							"to the following location: %s\n", filepath.c_str());
			throw(runtime_error(message));
		}
		while (!dataFile.eof()) {
			std::string temp = "";
			while (std::getline(dataFile, temp)) {
				tabidx = temp.find('\t');
				return_idx = temp.find('\r');
				if (!return_idx)
					return_idx = temp.length();
				if (tabidx >= 0 && tabidx < (int)temp.length()) {
					v1 = temp.substr(0, tabidx);
					v2 = temp.substr(tabidx + 1, return_idx - tabidx - 1);
					//std::cout << "row" << v1 << ", " << v2 << std::endl;
					mapping->insert( { v1, v2 });
				}
				temp = "";
			}
		}
		return mapping;

	}

	virtual ~Parser();
};

#endif /* PARSER_H_ */

