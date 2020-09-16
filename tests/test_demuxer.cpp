#define BOOST_TEST_MODULE boost_test_barcode
#include <boost/test/included/unit_test.hpp>
#include <boost/filesystem/path.hpp>

#include <unordered_map>
#include <string>
#include "../src/Parser.h"
#include "../src/Barcode.h"
#include "../src/PairedReader.h"
#include "../src/Demux.h"
#include "../src/helper.h"

// compile with:  g++ -std=c++11 test_demuxer.cpp ../src/helper.h ../src/Barcode.h ../src/Barcode.cpp ../src/Parser.h ../src/Parser.cpp ../src/PairedReader.h ../src/PairedReader.cpp ../src/ZipFastqReader.h ../src/ZipFastqReader.cpp -lboost_system -lboost_filesystem -lpthread -lz -ldl

using namespace std;

int I1_START = 10;

pair<int, int> I7_POS = { 0, 12 };
pair<int, int> I5_POS = { 12, 24 };
pair<int, int> I1_POS = { 24, 36 };

string Exe_path = utils::getExecutablePath().length() > 0 ? boost::filesystem::path(utils::getExecutablePath()).parent_path().string() : std::string();

vector<string>* demux_i7_i5_i1() {
	string res = Exe_path + PATH_SEP + string("resources") + PATH_SEP
			+ string("end_to_end");
	string read_1 = res + PATH_SEP + string("i7_i5_i1_read_1.fastq.gz");
	string read_2 = res + PATH_SEP + string("i7_i5_i1_read_2.fastq.gz");
	string csv = res + PATH_SEP + string("i7_i5_i1_sample_sheet.csv");
	vector<string> *paths = new vector<string>( { read_1, read_2, csv });
	return paths;
}

vector<string>* demux_i7_i1() {
	string res = Exe_path + PATH_SEP + string("resources") + PATH_SEP
			+ string("end_to_end");
	string read_1 = res + PATH_SEP + string("i7_i1_read_1.fastq.gz");
	string read_2 = res + PATH_SEP + string("i7_i1_read_2.fastq.gz");
	string csv = res + PATH_SEP + string("i7_i1_sample_sheet.csv");
	vector<string> *paths = new vector<string>( { read_1, read_2, csv });
	return paths;
}

vector<string>* demux_i5_i1() {
	string res = Exe_path + PATH_SEP + string("resources") + PATH_SEP
			+ string("end_to_end");
	string read_1 = res + PATH_SEP + string("i5_i1_read_1.fastq.gz");
	string read_2 = res + PATH_SEP + string("i5_i1_read_2.fastq.gz");
	string csv = res + PATH_SEP + string("i5_i1_sample_sheet.csv");
	vector<string> *paths = new vector<string>( { read_1, read_2, csv });
	return paths;
}

vector<string>* demux_i1() {
	// demulitplexing only on i1 one.
	// should not matter if the header contains i7 and i5 barcodes
	string res = Exe_path + PATH_SEP + string("resources") + PATH_SEP
			+ string("end_to_end");
	string read_1 = res + PATH_SEP + string("i7_i5_i1_read_1.fastq.gz");
	string read_2 = res + PATH_SEP + string("i7_i5_i1_read_2.fastq.gz");
	string csv = res + PATH_SEP + string("i1_sample_sheet.csv");
	vector<string> *paths = new vector<string>( { read_1, read_2, csv });
	return paths;
}

vector<string>* demux_i7_i5() {
	string res = Exe_path + PATH_SEP + string("resources") + PATH_SEP
			+ string("end_to_end");
	string read_1 = res + PATH_SEP + string("i7_i5_read_1.fastq.gz");
	string read_2 = res + PATH_SEP + string("i7_i5_read_2.fastq.gz");
	string csv = res + PATH_SEP + string("i7_i5_sample_sheet.csv");
	vector<string> *paths = new vector<string>( { read_1, read_2, csv });
	return paths;
}

void demux_loop(string read1, string read2, vector<Barcode*> &barcodes,
		int i1_start, int i1_end,
		void cb(string, std::pair<fq_read*, fq_read*>&)) {
	//i7, i5, i1 = barcodes
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
	PairedReader pr(read1, read2);
	std::vector<std::pair<fq_read*, fq_read*>> *pe_reads;
	bool reads_available = true;
	std::pair<fq_read*, fq_read*> processed_mates;
	while (reads_available) {
		pe_reads = pr.next_reads2(40000, 2);
		if (pe_reads == NULL || pe_reads->size() == 0)
			reads_available = false;
		//here we do the error correction and get obtain the i1 barcode if present
		for (auto it_mp = pe_reads->begin(); it_mp != pe_reads->end();
				it_mp++) {
			pair<fq_read*, fq_read*> mate_pair = *it_mp;
			string s_barcodes = process_mate_pair(mate_pair, i7_wanted,
					i5_wanted, i1_wanted, map_i7, map_i5, map_i1, i1_start,
					i1_end, processed_mates);
			cb(s_barcodes, processed_mates);
			delete processed_mates.first;
			delete processed_mates.second;
		}
		delete pe_reads;
	}
}

void test_demux_i7_i5_i1_callback(string corrected_bc,
		std::pair<fq_read*, fq_read*> &processed_mates) {
	int first_sep = corrected_bc.find_first_of("\n");
	int second_sep = corrected_bc.find_first_of("\n", first_sep + 1);
	int third_sep = corrected_bc.find_first_of("\n", second_sep + 1);
	string i7_corr = corrected_bc.substr(0, first_sep);
	string i5_corr = corrected_bc.substr(first_sep + 1,
			second_sep - first_sep - 1);
	string i1_corr = corrected_bc.substr(second_sep + 1,
			corrected_bc.length() - second_sep - 1);
	fq_read *m1 = processed_mates.first;
	fq_read *m2 = processed_mates.second;

	string expected_i7 = m1->Sequence.substr(I7_POS.first,
			I7_POS.second - I7_POS.first);
	string expected_i5 = m1->Sequence.substr(I5_POS.first,
			I5_POS.second - I5_POS.first);
	string expected_i1 = m1->Sequence.substr(I1_POS.first,
			I1_POS.second - I1_POS.first);

	BOOST_CHECK(i7_corr == expected_i7);
	//std::cout << i5_corr << " " << i5_corr << std::endl;
	BOOST_CHECK(i5_corr == expected_i5);
	BOOST_CHECK(i1_corr == expected_i1);

	BOOST_CHECK(m1->Sequence.length() == m1->QualityCode.length());
	BOOST_CHECK(m2->Sequence.length() == m2->QualityCode.length());
}
/**
 * Testing expected vs calcualted barcodes from paired-end fastq files

 * reads in read1 fastq file look as follows:
 <read_name>:<correctable_i7>+<correctable_i5>
 <correct_i7><correct_i5><correct_i1>
 +
 E{36}

 * reads in read2 fastq file look as follows:
 <read_name>:<correctable_i7>+<correctable_i5>
 <rand_10nt><correctable_i1>
 +
 E{22}
 */
BOOST_AUTO_TEST_CASE( test_demux_i7_i5_i1 ) {
	vector<string> *paths = demux_i7_i5_i1();

	string read1, read2, csv;
	read1 = paths->at(0);
	read2 = paths->at(1);
	csv = paths->at(2);
	Parser pe;
	vector<Barcode*> barcodes;
	std::cout << csv << std::endl;
	boost::filesystem::path p(utils::getExecutablePath()); //.parent_path();
	unordered_map<string, string> *barcode_sample_map = pe.parse_sample_sheet(
			csv, false, barcodes, p.string());

	Barcode *i1 = barcodes[2];

	int i1_end = I1_START + i1->length;

	demux_loop(read1, read2, barcodes, I1_START, i1_end,
			test_demux_i7_i5_i1_callback);
}

void test_demux_i7_i1_callback(string corrected_bc,
		std::pair<fq_read*, fq_read*> &processed_mates) {
	int first_sep = corrected_bc.find_first_of("\n");
	int second_sep = corrected_bc.find_first_of("\n", first_sep + 1);
	int third_sep = corrected_bc.find_first_of("\n", second_sep + 1);
	string i7_corr = corrected_bc.substr(0, first_sep);
	string i5_corr = corrected_bc.substr(first_sep + 1,
			second_sep - first_sep - 1);
	string i1_corr = corrected_bc.substr(second_sep + 1,
			corrected_bc.length() - second_sep - 1);
	fq_read *m1 = processed_mates.first;
	fq_read *m2 = processed_mates.second;

	string expected_i7 = m1->Sequence.substr(I7_POS.first,
			I7_POS.second - I7_POS.first);
	string expected_i5 = "";
	string expected_i1 = m1->Sequence.substr(I1_POS.first,
			I1_POS.second - I1_POS.first);

	BOOST_CHECK(i7_corr == expected_i7);
	//std::cout << i5_corr << " " << i5_corr << std::endl;
	BOOST_CHECK(i5_corr == expected_i5);
	BOOST_CHECK(i1_corr == expected_i1);

	BOOST_CHECK(m1->Sequence.length() == m1->QualityCode.length());
	BOOST_CHECK(m2->Sequence.length() == m2->QualityCode.length());
}

/**
 * Testing expected vs calcualted barcodes from paired-end fastq files

 * reads in read1 fastq file look as follows:
 <read_name>:<correctable_i7>+<correctable_i5>
 <correct_i7><correct_i5><correct_i1>
 +
 E{36}

 * reads in read2 fastq file look as follows:
 <read_name>:<correctable_i7>+<correctable_i5>
 <rand_10nt><correctable_i1>
 +
 E{22}
 */
BOOST_AUTO_TEST_CASE(test_demux_i7_i1) {
	vector<string> *paths = demux_i7_i1();

	string read1, read2, csv;
	read1 = paths->at(0);
	read2 = paths->at(1);
	csv = paths->at(2);
	Parser pe;
	vector<Barcode*> barcodes;
	std::cout << csv << std::endl;
	boost::filesystem::path p(utils::getExecutablePath()); //.parent_path();
	unordered_map<string, string> *barcode_sample_map = pe.parse_sample_sheet(
			csv, false, barcodes, p.string());
	Barcode *i1 = barcodes[2];

	int i1_end = I1_START + i1->length;

	demux_loop(read1, read2, barcodes, I1_START, i1_end,
			test_demux_i7_i1_callback);
}

void test_demux_i7_i5_callback(string corrected_bc,
		std::pair<fq_read*, fq_read*> &processed_mates) {
	int first_sep = corrected_bc.find_first_of("\n");
	int second_sep = corrected_bc.find_first_of("\n", first_sep + 1);
	int third_sep = corrected_bc.find_first_of("\n", second_sep + 1);
	string i7_corr = corrected_bc.substr(0, first_sep);
	string i5_corr = corrected_bc.substr(first_sep + 1,
			second_sep - first_sep - 1);
	string i1_corr = corrected_bc.substr(second_sep + 1,
			corrected_bc.length() - second_sep - 1);
	fq_read *m1 = processed_mates.first;
	fq_read *m2 = processed_mates.second;

	string expected_i7 = m1->Sequence.substr(I7_POS.first,
			I7_POS.second - I7_POS.first);
	string expected_i5 = m1->Sequence.substr(I5_POS.first,
			I5_POS.second - I5_POS.first);
	string expected_i1 = "";

	BOOST_CHECK(i7_corr == expected_i7);
	//std::cout << i5_corr << " " << i5_corr << std::endl;
	BOOST_CHECK(i5_corr == expected_i5);
	BOOST_CHECK(i1_corr == expected_i1);

	BOOST_CHECK(m1->Sequence.length() == m1->QualityCode.length());
	BOOST_CHECK(m2->Sequence.length() == m2->QualityCode.length());
}

BOOST_AUTO_TEST_CASE(test_demux_i7_i5) {
	vector<string> *paths = demux_i7_i5();

	string read1, read2, csv;
	read1 = paths->at(0);
	read2 = paths->at(1);
	csv = paths->at(2);
	Parser pe;
	vector<Barcode*> barcodes;
	std::cout << csv << std::endl;
	boost::filesystem::path p(utils::getExecutablePath()); //.parent_path();
	unordered_map<string, string> *barcode_sample_map = pe.parse_sample_sheet(
			csv, false, barcodes, p.string());
	Barcode *i1 = barcodes[2];

	int i1_end = I1_START + i1->length;

	demux_loop(read1, read2, barcodes, I1_START, i1_end,
			test_demux_i7_i5_callback);
}

void test_demux_i1_callback(string corrected_bc,
		std::pair<fq_read*, fq_read*> &processed_mates) {
	int first_sep = corrected_bc.find_first_of("\n");
	int second_sep = corrected_bc.find_first_of("\n", first_sep + 1);
	int third_sep = corrected_bc.find_first_of("\n", second_sep + 1);
	string i7_corr = corrected_bc.substr(0, first_sep);
	string i5_corr = corrected_bc.substr(first_sep + 1,
			second_sep - first_sep - 1);
	string i1_corr = corrected_bc.substr(second_sep + 1,
			corrected_bc.length() - second_sep - 1);
	fq_read *m1 = processed_mates.first;
	fq_read *m2 = processed_mates.second;

	string expected_i7 = "";
	string expected_i5 = "";
	string expected_i1 = m1->Sequence.substr(I1_POS.first,
			I1_POS.second - I1_POS.first);

	BOOST_CHECK(i7_corr == expected_i7);
	BOOST_CHECK(i5_corr == expected_i5);
	BOOST_CHECK(i1_corr == expected_i1);

	BOOST_CHECK(m1->Sequence.length() == m1->QualityCode.length());
	BOOST_CHECK(m2->Sequence.length() == m2->QualityCode.length());
}

BOOST_AUTO_TEST_CASE(test_demux_i1) {
	vector<string> *paths = demux_i1();

	string read1, read2, csv;
	read1 = paths->at(0);
	read2 = paths->at(1);
	csv = paths->at(2);
	Parser pe;
	vector<Barcode*> barcodes;
	std::cout << csv << std::endl;
	boost::filesystem::path p(utils::getExecutablePath()); //.parent_path();
	unordered_map<string, string> *barcode_sample_map = pe.parse_sample_sheet(
			csv, false, barcodes, p.string());
	Barcode *i1 = barcodes[2];

	int i1_end = I1_START + i1->length;

	demux_loop(read1, read2, barcodes, I1_START, i1_end,
			test_demux_i1_callback);
}

BOOST_AUTO_TEST_CASE(test_demux_paired_end) {
	int expected_reads = 100;

	vector<string> *paths = demux_i7_i5_i1();
	string read1, read2, csv;
	read1 = paths->at(0);
	read2 = paths->at(1);
	csv = paths->at(2);
	Parser pe;
	vector<Barcode*> barcodes;
	std::cout << csv << std::endl;
	boost::filesystem::path p(utils::getExecutablePath()); //.parent_path();
	unordered_map<string, string> *barcode_sample_map = pe.parse_sample_sheet(
			csv, false, barcodes, p.string());

	int chunk_size = 1000;
	int reading_threads = 2;
	int writing_threads = 1;
	int processing_threads = 1;

	string tmp_path = p.parent_path().string() + PATH_SEP + string("test_output");
	if(!utils::folder_exists(tmp_path))
		utils::mkdir(tmp_path);
	demux_paired_end(barcode_sample_map, barcodes, read1, read2, I1_START,
			tmp_path, pe, chunk_size, reading_threads, writing_threads,
			processing_threads);

	//string stats_file = tmp_path + PATH_SEP + string("demultipexing_stats.tsv");
	//ifstream stats_file_stream(stats_file.c_str());
	std::unordered_map<string, string> *tsv_lines =
			Parser::get_map_from_resource(tmp_path, "demultipexing_stats.tsv");
	string read_count;
	for (auto it = tsv_lines->begin(); it != tsv_lines->end(); it++) {
		// test all sample lines (without header line).
		if (it->first.rfind("sample_", 0) == 0) {
			read_count = it->second;
			stringstream ss(read_count);
			int n_reads;
			ss >> n_reads; //strtol(row[1].c_str(), &rest, 10);
			BOOST_CHECK(n_reads == expected_reads);
		}
	}
	delete tsv_lines;
	utils::rmdir(tmp_path);
}

// TODO write test for mixed cases

