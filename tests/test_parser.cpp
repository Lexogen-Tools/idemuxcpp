#define BOOST_TEST_MODULE boost_test_barcode
#include <boost/test/included/unit_test.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/filesystem/path.hpp>

#include <iostream>
#include <unordered_map>
#include <string>
#include "../src/Parser.h"
#include "../src/Barcode.h"
#include "../src/helper.h"
#include "../src/PairedReader.h"
#include <stdexcept>

// compile with:  g++ -std=c++11 test_parser.cpp ../src/helper.h ../src/Barcode.h ../src/Barcode.cpp ../src/Parser.h ../src/Parser.cpp ../src/PairedReader.h ../src/PairedReader.cpp ../src/ZipFastqReader.h ../src/ZipFastqReader.cpp -lboost_system -lboost_filesystem -lpthread -lz -ldl

using namespace std;

string RESOURCE = "resources";

string Exe_path = utils::getExecutablePath();
string Exe_dir = Exe_path.length() > 0 ? boost::filesystem::path(Exe_path).parent_path().string() : std::string();

string CSVS_TO_FAIL = Exe_dir + PATH_SEP + RESOURCE + PATH_SEP + string("sample_sheet/fail/*.csv");
string CSVS_TO_PASS = Exe_dir + PATH_SEP + RESOURCE + PATH_SEP + string("sample_sheet/pass/*.csv");

vector<Barcode*> BARCODES_TO_TEST() {
	string A6 = string(6, 'A');
	string C6 = string(6, 'C');
	string A12 = string(12, 'A');
	string C12 = string(12, 'C');
	unordered_map<string, vector<string>> BARCODES_6NT;
	BARCODES_6NT[A6] = { "test0" };
	BARCODES_6NT[C6] = { "test1" };

	unordered_map<string, vector<string>> BARCODES_NON_LEX;
	BARCODES_NON_LEX[A12] = { "test0" };
	BARCODES_NON_LEX[C12] = { "test1" };

	unordered_map<string, vector<string>> BARCODES_NONE;
	BARCODES_NONE[""] = { "test0", "test1" };

	vector<Barcode*> bc;
	bc.push_back(new Barcode("i7", BARCODES_6NT));
	bc.push_back(new Barcode("i7", BARCODES_NON_LEX));
	bc.push_back(new Barcode("i7", BARCODES_NONE));
	return bc;
}

vector<pair<string, string>> FOR_RC_PASS = { { "AAAACATCGTTN", "NAACGATGTTTT" },
		{ "", "" }, { "", "" } };
vector<string> FOR_RC_FAIL = { "AAAACATCGTTNXXXXX", to_string(10) };

string FQ_RES = Exe_dir + PATH_SEP + RESOURCE + PATH_SEP + string("fastq");
int INDEX_LENGTH = 12;
int I1_START = 10;

struct fq_test_data {
	string condition;
	string fq_gz_1;
	string fq_gz_2;
	bool has_i7;
	bool has_i5;
	bool has_i1;
	int i7_length;
	int i5_length;
	int i1_start;
	int i1_end;
};

vector<fq_test_data> FQ_TO_PASS = { { "i7(12)nt, i5(12)nt, i1(12)nt", FQ_RES
		+ PATH_SEP + string("i7_i5_i1_read_1.fastq.gz"), FQ_RES + PATH_SEP
		+ string("i7_i5_i1_read_2.fastq.gz"), true, true, true, INDEX_LENGTH,
		INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH }, {
		"i7(12)nt, i5(0nt), i1(12)nt", FQ_RES + PATH_SEP
				+ string("i7_i1_read_1.fastq.gz"), FQ_RES + PATH_SEP
				+ string("i7_i1_read_2.fastq.gz"), true, false, true,
		INDEX_LENGTH, 0, I1_START, I1_START + INDEX_LENGTH }, {
		"i7(0)nt, i5(12nt), i1(12)nt", FQ_RES + PATH_SEP
				+ string("i5_i1_read_1.fastq.gz"), FQ_RES + PATH_SEP
				+ string("i5_i1_read_2.fastq.gz"), false, true, true, 0,
		INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH }, {
		"i7(10)nt, i5(12nt), i1(12)nt", FQ_RES + PATH_SEP
				+ string("i7-2_i5_i1_read_1.fastq.gz"), FQ_RES + PATH_SEP
				+ string("i7-2_i5_i1_read_2.fastq.gz"), true, true, true,
		INDEX_LENGTH - 2, INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH }, };

vector<fq_test_data> FQ_TO_FAIL = { { "too long i7", FQ_RES + PATH_SEP
		+ string("i7+2_i5_i1_read_1.fastq.gz"), FQ_RES + PATH_SEP
		+ string("i7_i5_i1_read_2.fastq.gz"), true, true, true, INDEX_LENGTH,
		INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH }, { "too long i5",
		FQ_RES + PATH_SEP + string("i7_i5+2_i1_read_1.fastq.gz"), FQ_RES
				+ PATH_SEP + string("i7_i5+2_i1_read_2.fastq.gz"), true, true,
		true, INDEX_LENGTH, INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH }, {
		"too short i1 sequence", FQ_RES + PATH_SEP
				+ string("i7_i5_noi1_read_1.fastq.gz"), FQ_RES + PATH_SEP
				+ string("i7_i5_noi1_read_2.fastq.gz"), true, true, true,
		INDEX_LENGTH, INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH }, {
		"different barcode headers", FQ_RES + PATH_SEP
				+ string("i7+2_i5_i1_read_1.fastq.gz"), FQ_RES + PATH_SEP
				+ string("i7_i5_i1_read_2.fastq.gz"), true, true, true,
		INDEX_LENGTH, INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH }, {
		"too short i7", FQ_RES + PATH_SEP
				+ string("i7-2_i5_i1_read_1.fastq.gz"), FQ_RES + PATH_SEP
				+ string("i7-2_i5_i1_read_2.fastq.gz"), true, true, true,
		INDEX_LENGTH, INDEX_LENGTH, I1_START, I1_START + INDEX_LENGTH } };

vector<string> iterate_csv_files(boost::filesystem::path p) {
	vector<string> csv_files;
	if (boost::filesystem::is_directory(p)) {
		//std::cout << p.string() << " is a directory containing:\n";
		for (auto &entry : boost::make_iterator_range(
				boost::filesystem::directory_iterator(p), { })) {

			string path = entry.path().string();

			if (path.rfind(".csv", path.length() - 4) != string::npos) {
				//std::cout << path << "\n";
				csv_files.push_back(path);
			}
		}
	}
	return csv_files;
}

boost::filesystem::path tmp_path(CSVS_TO_FAIL);
vector<string> csv_files_to_fail = iterate_csv_files(tmp_path.parent_path());

BOOST_AUTO_TEST_CASE( test_parse_sample_sheet_to_fail ) {
	Parser pe;
	vector<Barcode*> barcodes;
	boost::filesystem::path p(Exe_path);
	for (auto it = csv_files_to_fail.begin(); it != csv_files_to_fail.end();
			it++) {
		string csv_file = *it;
		//std::cout << csv_file << std::endl;
		unordered_map<string, string> *rows = NULL;
		BOOST_CHECK_THROW(
				rows = pe.parse_sample_sheet(csv_file, false, barcodes,
						p.string()), std::exception);
		if (rows)
			delete rows;
	}
	for (auto it = barcodes.begin(); it != barcodes.end(); it++)
		delete *it;
}

boost::filesystem::path tmp_path2(CSVS_TO_PASS);
vector<string> csv_files_to_pass = iterate_csv_files(tmp_path2.parent_path());

BOOST_AUTO_TEST_CASE( test_parse_sample_sheet_to_pass ) {
	Parser pe;
	vector<Barcode*> barcodes;
	boost::filesystem::path p(Exe_path);
	for (auto it = csv_files_to_pass.begin(); it != csv_files_to_pass.end();
			it++) {
		string csv_file = *it;
		//std::cout << csv_file << std::endl;
		BOOST_CHECK_NO_THROW(
				pe.parse_sample_sheet(csv_file, false, barcodes, p.string()));
	}
}

BOOST_AUTO_TEST_CASE( test_load_correction_maps ) {
	Parser pe;
	vector<Barcode*> bctest = BARCODES_TO_TEST();
	for (auto it = bctest.begin(); it != bctest.end(); it++) {
		Barcode *tmp_bc = *it;
		boost::filesystem::path p(Exe_path);
		tmp_bc->load_correction_map(p.string());
		if (!tmp_bc->correction_map) {
			throw(std::runtime_error("ERROR: could not load correction map!\n"));
		}
		for (auto itk = tmp_bc->correction_map->begin();
				itk != tmp_bc->correction_map->end(); itk++) {
			BOOST_CHECK(itk->first.compare(itk->second) == 0);
		}

		delete tmp_bc;
	}
}

BOOST_AUTO_TEST_CASE( test_reverse_complement_pass ) {
	Parser p;
	for (auto it = FOR_RC_PASS.begin(); it != FOR_RC_PASS.end(); it++) {
		string test_in = it->first;
		string expected = it->second;
		string res = p.reverse_complement(test_in);
		BOOST_CHECK(res.compare(expected) == 0);
	}
}

BOOST_AUTO_TEST_CASE( test_reverse_complement_fail ) {
	Parser p;
	for (auto it = FOR_RC_FAIL.begin(); it != FOR_RC_FAIL.end(); it++) {
		string test_in = *it;
		BOOST_CHECK_THROW(p.reverse_complement(test_in), std::exception);
	}
}

BOOST_AUTO_TEST_CASE( test_peek_into_fastq_files_fq_pass ) {
	Parser p;
	for (auto it = FQ_TO_PASS.begin(); it != FQ_TO_PASS.end(); it++) {
		BOOST_CHECK_NO_THROW(
				p.peek_into_fastq_files(it->fq_gz_1, it->fq_gz_2, it->has_i7,
						it->has_i5, it->has_i1, it->i7_length, it->i5_length,
						it->i1_start, it->i1_end));
	}
}

BOOST_AUTO_TEST_CASE( test_peek_into_fastq_files_fq_fail ) {
	Parser p;
	for (auto it = FQ_TO_FAIL.begin(); it != FQ_TO_FAIL.end(); it++) {
		BOOST_CHECK_THROW(
				p.peek_into_fastq_files(it->fq_gz_1, it->fq_gz_2, it->has_i7,
						it->has_i5, it->has_i1, it->i7_length, it->i5_length,
						it->i1_start, it->i1_end), std::exception);
	}
}

bool is_equal_read(fq_read &r1, fq_read r2) {
	bool is_equal = r1.Seq_ID.compare(r2.Seq_ID) == 0
			&& r1.Sequence.compare(r2.Sequence) == 0
			&& r1.Plus_ID.compare(r2.Plus_ID) == 0
			&& r1.QualityCode.compare(r2.QualityCode) == 0;
	return is_equal;
}

BOOST_AUTO_TEST_CASE( test_fq_gz_parser ) {
	fq_read test_r_1, test_r_2;
	test_r_1.Seq_ID = "@read_name:AANACATGCGTT+CGCCACTGAGTT";
	test_r_1.Sequence = "AAAACATGCGTTCCCCACTGAGTTAAAACATGCGTT";
	test_r_1.Plus_ID = "+";
	test_r_1.QualityCode = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE";
	test_r_2.Seq_ID = "@read_name:AANACATGCGTT+CGCCACTGAGTT";
	test_r_2.Sequence = "TTTTAGCATGAAAACATGCGTT";
	test_r_2.Plus_ID = "+";
	test_r_2.QualityCode = "EEEEEEEEEEEEEEEEEEEEEE";

	string fq_1 = FQ_RES + PATH_SEP + string("i7_i5_i1_read_1.fastq.gz");
	string fq_2 = FQ_RES + PATH_SEP + string("i7_i5_i1_read_2.fastq.gz");
	PairedReader get_pe_fastq(fq_1, fq_2);
	std::vector<std::pair<fq_read*, fq_read*>> *pe_reads =
			get_pe_fastq.next_reads(100);
	for (auto it = pe_reads->begin(); it != pe_reads->end(); it++) {
		BOOST_CHECK(is_equal_read(*it->first, test_r_1) == true);
		BOOST_CHECK(is_equal_read(*it->second, test_r_2) == true);
		delete it->first;
		delete it->second;
	}
	delete pe_reads;
}
