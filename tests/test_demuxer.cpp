#define BOOST_TEST_MODULE boost_test_demuxer
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

int I1_START = 10; // zero based.

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

vector<string>* demux_i7() {
	string res = Exe_path + PATH_SEP + string("resources") + PATH_SEP
			+ string("end_to_end");
	string read_1 = res + PATH_SEP + string("test_12_10_6_r1.fastq.gz");
	string read_2 = res + PATH_SEP + string("test_12_10_6_r2.fastq.gz");
	string csv = res + PATH_SEP + string("i7_different_lengths.csv");
	vector<string> *paths = new vector<string>( { read_1, read_2, csv });
	return paths;
}

vector<string>* demux_i7_i1_different_lengths(){
	string res = Exe_path + PATH_SEP + string("resources") + PATH_SEP
			+ string("end_to_end");
	string read_1 = res + PATH_SEP + string("test_12_10_6_r1.fastq.gz");
	string read_2 = res + PATH_SEP + string("test_12_10_6_r2.fastq.gz");
	string csv = res + PATH_SEP + string("i7_i1_different_lengths.csv");
	vector<string> *paths = new vector<string>( { read_1, read_2, csv });
	return paths;
}

vector<string>* demux_i7_i1_different_lengths_se(){
	string res = Exe_path + PATH_SEP + string("resources") + PATH_SEP
			+ string("end_to_end");
	string read_1 = res + PATH_SEP + string("test_12_10_6_r1.fastq.gz");
	string read_2 = res + PATH_SEP + string("test_12_10_6_r2.fastq.gz");
	string csv = res + PATH_SEP + string("i7_i1_different_lengths_se.csv");
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
		unordered_map<string, i1_info> &i7_i5_i1_info_map,
		void cb(string, std::pair<fq_read*, fq_read*>&)) {
	//i7, i5, i1 = barcodes
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
	PairedReader pr(read1, read2, 0);
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
					i5_wanted, i1_wanted, map_i7, map_i5, map_i1, i7_i5_i1_info_map, processed_mates, NULL);
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
	unordered_map<string, i1_info> i7_i5_i1_info_map;
	std::cout << csv << std::endl;
	boost::filesystem::path p(utils::getExecutablePath()); //.parent_path();
	unordered_map<string, string> *barcode_sample_map = pe.parse_sample_sheet(
			csv, false, barcodes,i7_i5_i1_info_map, p.string(), "", false, 2, I1_START);

	Barcode *i1 = barcodes[2];

	//int i1_end = I1_START + i1->length;

	demux_loop(read1, read2, barcodes, i7_i5_i1_info_map,
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
	unordered_map<string, i1_info> i7_i5_i1_info_map;
	std::cout << csv << std::endl;
	boost::filesystem::path p(utils::getExecutablePath()); //.parent_path();
	unordered_map<string, string> *barcode_sample_map = pe.parse_sample_sheet(
			csv, false, barcodes, i7_i5_i1_info_map, p.string(), "", false, 2, I1_START);
	Barcode *i1 = barcodes[2];

	//int i1_end = I1_START + i1->length;

	demux_loop(read1, read2, barcodes, i7_i5_i1_info_map,
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
	unordered_map<string, i1_info> i7_i5_i1_info_map;
	std::cout << csv << std::endl;
	boost::filesystem::path p(utils::getExecutablePath()); //.parent_path();
	unordered_map<string, string> *barcode_sample_map = pe.parse_sample_sheet(
			csv, false, barcodes, i7_i5_i1_info_map, p.string(), "", false, 2, I1_START);
	Barcode *i1 = barcodes[2];

	//int i1_end = I1_START + i1->length;

	demux_loop(read1, read2, barcodes, i7_i5_i1_info_map,
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
	unordered_map<string, i1_info> i7_i5_i1_info_map;
	std::cout << csv << std::endl;
	boost::filesystem::path p(utils::getExecutablePath()); //.parent_path();
	unordered_map<string, string> *barcode_sample_map = pe.parse_sample_sheet(
			csv, false, barcodes, i7_i5_i1_info_map, p.string(), "", false, 2, I1_START);
	Barcode *i1 = barcodes[2];

	//int i1_end = I1_START + i1->length;

	demux_loop(read1, read2, barcodes, i7_i5_i1_info_map,
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
	unordered_map<string, i1_info> i7_i5_i1_info_map;
	std::cout << csv << std::endl;
	boost::filesystem::path p(utils::getExecutablePath()); //.parent_path();
	unordered_map<string, string> *barcode_sample_map = pe.parse_sample_sheet(
			csv, false, barcodes, i7_i5_i1_info_map, p.string(), "", false, 2, I1_START);

	int chunk_size = 1000;
	int reading_threads = 2;
	int writing_threads = 1;
	int processing_threads = 1;

	string tmp_path = p.parent_path().string() + PATH_SEP + string("test_output");
	if(!utils::folder_exists(tmp_path))
		utils::mkdir(tmp_path);
	bool skip_check = false;
	bool restrict_barcode_length = false;
	double size_writer_buffer_gb = 1;
	double queue_buffer_gb = 0;
	bool verbose = true;
	demux_paired_end(barcode_sample_map, barcodes, read1, read2, i7_i5_i1_info_map,
			tmp_path, pe, chunk_size, reading_threads, writing_threads,
			processing_threads, tmp_path + PATH_SEP + string("count_corrections.tsv"), skip_check, restrict_barcode_length, size_writer_buffer_gb, queue_buffer_gb, verbose);

	//string stats_file = tmp_path + PATH_SEP + string("demultipexing_stats.tsv");
	//ifstream stats_file_stream(stats_file.c_str());
	std::unordered_map<string, string> *tsv_lines =
			Parser::get_map_from_resource(tmp_path, "demultipexing_stats.tsv");
	string read_count;
	for (auto it = tsv_lines->begin(); it != tsv_lines->end(); it++) {
		// test all sample lines (without header line).
		if (it->first.rfind("sample_", 0) == 0 && it->first.compare("sample_name") != 0) {
			read_count = it->second;
			stringstream ss(read_count);
			int n_reads;
			ss >> n_reads; //strtol(row[1].c_str(), &rest, 10);
			fprintf(stderr,"reads %d 100 %s\n", n_reads, it->first.c_str());
			BOOST_CHECK(n_reads == expected_reads);
		}
	}
	delete tsv_lines;
	utils::rmdir(tmp_path);
}

void test_demux_reads(vector<string> *paths_r1_r2_csv){
	int expected_reads = 1;
	string read1, read2, csv;
	read1 = paths_r1_r2_csv->at(0);
	read2 = paths_r1_r2_csv->at(1);
	csv = paths_r1_r2_csv->at(2);
	Parser pe;
	vector<Barcode*> barcodes;
	unordered_map<string, i1_info> i7_i5_i1_info_map;
	std::cout << csv << std::endl;
	boost::filesystem::path p(utils::getExecutablePath()); //.parent_path();
	unordered_map<string, string> *barcode_sample_map = pe.parse_sample_sheet(
			csv, false, barcodes, i7_i5_i1_info_map, p.string(), "", false, 2, I1_START);

	int chunk_size = 1000;
	int reading_threads = 2;
	int writing_threads = 1;
	int processing_threads = 1;

	string tmp_path = p.parent_path().string() + PATH_SEP + string("test_output");
	if(!utils::folder_exists(tmp_path))
		utils::mkdir(tmp_path);
	bool skip_check = false;
	bool restrict_barcode_length = false;
	double size_writer_buffer_gb = 1;
	double queue_buffer_gb = 0;
	bool verbose = true;
	demux_paired_end(barcode_sample_map, barcodes, read1, read2, i7_i5_i1_info_map,
		tmp_path, pe, chunk_size, reading_threads, writing_threads,
		processing_threads, tmp_path + PATH_SEP + string("count_corrections.tsv"), skip_check, restrict_barcode_length, size_writer_buffer_gb, queue_buffer_gb, verbose);

	//string stats_file = tmp_path + PATH_SEP + string("demultipexing_stats.tsv");
	//ifstream stats_file_stream(stats_file.c_str());
	std::unordered_map<string, string> *tsv_lines =
			Parser::get_map_from_resource(tmp_path, "demultipexing_stats.tsv");
	string read_count;
	for (auto it = tsv_lines->begin(); it != tsv_lines->end(); it++) {
		// test all sample lines (without header line).
		if (it->first.rfind("sample_", 0) == 0 && it->first.compare("sample_name") != 0) {
			read_count = it->second;
			stringstream ss(read_count);
			int n_reads;
			ss >> n_reads; //strtol(row[1].c_str(), &rest, 10);
			fprintf(stderr,"reads %d 1 %s\n", n_reads, it->first.c_str());
			BOOST_CHECK(n_reads == expected_reads);
		}
	}
	delete tsv_lines;
	utils::rmdir(tmp_path);
}

BOOST_AUTO_TEST_CASE(test_demux_paired_end_i7) {
	vector<string> *paths = demux_i7();
	test_demux_reads(paths);
}

BOOST_AUTO_TEST_CASE(test_demux_paired_end_i7_i1) {
	vector<string> *paths = demux_i7_i1_different_lengths();
	test_demux_reads(paths);
}

BOOST_AUTO_TEST_CASE(test_demux_single_end_i7_i1) {
	vector<string> *paths = demux_i7_i1_different_lengths_se();
	test_demux_reads(paths);
}


bool is_equal_read(fq_read &r1, fq_read r2) {
	bool is_equal = r1.Seq_ID.compare(r2.Seq_ID) == 0;
	if(!is_equal)
		printf("Error: different id %s %s\n", r1.Seq_ID.c_str(), r2.Seq_ID.c_str());
	is_equal = r1.Sequence.compare(r2.Sequence) == 0 & is_equal;
	if(!is_equal)
			printf("Error: different sequence %s %s\n", r1.Sequence.c_str(), r2.Sequence.c_str());
	is_equal = r1.Plus_ID.compare(r2.Plus_ID) == 0 & is_equal;
	if(!is_equal)
			printf("Error: different plus id %s %s\n", r1.Plus_ID.c_str(), r2.Plus_ID.c_str());
	is_equal = r1.QualityCode.compare(r2.QualityCode) == 0 & is_equal;
	if(!is_equal)
			printf("Error: different quality code %s %s\n", r1.QualityCode.c_str(), r2.QualityCode.c_str());
	return is_equal;
}

BOOST_AUTO_TEST_CASE( test_demux_i7_i1_1_read ) {
	boost::filesystem::path p(utils::getExecutablePath()); //.parent_path();
	string tmp_path = p.parent_path().string() + PATH_SEP + string("test_output");
	if(!utils::folder_exists(tmp_path))
		utils::mkdir(tmp_path);

	string samplesheet = string("sample_name,i7,i5,i1\n") + string("test_sample,TCGAGTCC,,ATTCAACC");
	string tmp_path_samplesheet = tmp_path + PATH_SEP + string("samplesheet.csv");
	ofstream fs(tmp_path_samplesheet, std::ofstream::out);
	fs << samplesheet;
	fs.close();

	fq_read test_r_1, test_r_2;
	test_r_1.Seq_ID = "@NB502007:425:HGLGJBGXG:1:11101:11059:2649 1:N:0:TCGAGTCC";
	test_r_1.Sequence = "CCTGCTTCCTCCAACCCGACATGTGTACCTCAGCTTTTTCCCTCACTTGCATCAATAAAGCTTCTG";
	test_r_1.Plus_ID = "+";
	test_r_1.QualityCode = "A/A6AEEEEEEEA/EEAEEEEEEEEEE<EEEEEEEAEEEEEEEEAEEEAEEEEAEE/EEEEEEEE/";

	test_r_2.Seq_ID = "@NB502007:425:HGLGJBGXG:1:11101:11059:2649 2:N:0:TCGAGTCC";
	test_r_2.Sequence = "CTCTAGCTATAATCAACC";
	test_r_2.Plus_ID = "+";
	test_r_2.QualityCode = "6AAAA//E//E/A6E<AE";

	string tmp_path_r1 = tmp_path + PATH_SEP + string("read1.fastq");
	ofstream fsr1(tmp_path_r1, std::ofstream::out);
	fsr1 << test_r_1.to_string();
	fsr1.close();

	string tmp_path_r2 = tmp_path + PATH_SEP + string("read2.fastq");
	ofstream fsr2(tmp_path_r2, std::ofstream::out);
	fsr2 << test_r_2.to_string();
	fsr2.close();

	Parser pe;
	vector<Barcode*> barcodes;
	unordered_map<string, i1_info> i7_i5_i1_info_map;
	unordered_map<string, string> *barcode_sample_map = pe.parse_sample_sheet(
			tmp_path_samplesheet, false, barcodes, i7_i5_i1_info_map, p.string(), "", false, 2, I1_START);

	int chunk_size = 1000;
	int reading_threads = 2;
	int writing_threads = 1;
	int processing_threads = 1;

	// test too small and large buffer.
	for(int i = 0; i < 2; i++){
		bool skip_check = false;
		bool restrict_barcode_length = false;
		double size_writer_buffer_gb = i;
		double queue_buffer_gb = 0;
		bool verbose = true;
		demux_paired_end(barcode_sample_map, barcodes, tmp_path_r1, tmp_path_r2, i7_i5_i1_info_map,
					tmp_path, pe, chunk_size, reading_threads, writing_threads,
					processing_threads, tmp_path + PATH_SEP + string("count_corrections.tsv"), skip_check, restrict_barcode_length, size_writer_buffer_gb, queue_buffer_gb, verbose);

		string tmp_r1 = tmp_path + PATH_SEP + string("test_sample_R1.fastq.gz");
		string tmp_r2 = tmp_path + PATH_SEP + string("test_sample_R2.fastq.gz");
		BOOST_CHECK(boost::filesystem::exists(tmp_r1));
		BOOST_CHECK(boost::filesystem::exists(tmp_r1));

		fq_read exp_r_1, exp_r_2;
		exp_r_1.Seq_ID = "@NB502007:425:HGLGJBGXG:1:11101:11059:2649 1:N:0:TCGAGTCC+AATCAACC";
		exp_r_1.Sequence = "CCTGCTTCCTCCAACCCGACATGTGTACCTCAGCTTTTTCCCTCACTTGCATCAATAAAGCTTCTG";
		exp_r_1.Plus_ID = "+";
		exp_r_1.QualityCode = "A/A6AEEEEEEEA/EEAEEEEEEEEEE<EEEEEEEAEEEEEEEEAEEEAEEEEAEE/EEEEEEEE/";

		exp_r_2.Seq_ID = "@NB502007:425:HGLGJBGXG:1:11101:11059:2649 2:N:0:TCGAGTCC+AATCAACC";
		exp_r_2.Sequence = "CTCTAGCTAT";
		exp_r_2.Plus_ID = "+";
		exp_r_2.QualityCode = "6AAAA//E//";

		PairedReader get_pe_fastq(tmp_r1, tmp_r2, 0);
		std::vector<std::pair<fq_read*, fq_read*>> *pe_reads = get_pe_fastq.next_reads(100);
			size_t n_r = pe_reads->size();
			printf("nr %zu %s %s", n_r, tmp_r1.c_str(), tmp_r2.c_str());
		BOOST_CHECK(n_r == 1);
		for (auto it = pe_reads->begin(); it != pe_reads->end(); it++) {
			BOOST_CHECK(is_equal_read(*it->first, exp_r_1) == true);
			BOOST_CHECK(is_equal_read(*it->second, exp_r_2) == true);
			delete it->first;
			delete it->second;
		}
	}
}



// TODO write test for mixed cases

