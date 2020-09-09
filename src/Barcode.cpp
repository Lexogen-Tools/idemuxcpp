#include <vector>
#include <string>
#include <unordered_map>
#include "Barcode.h"
#include <iostream>
#include <istream>
#include <fstream>
#include "helper.h"
#include <stdexcept>
#include <algorithm>
#include "Parser.h"
using namespace std;

Barcode::Barcode(string barcode_type,  unordered_map<string,std::vector<string>> &ix_barcodes, bool reverse_complement) : length(0),Name(barcode_type),
		allowed_lengths({6, 8, 10, 12}),
    lengths_96_barcodes({8, 10, 12}),
    lengths_384_barcodes({10, 12}), correction_map(NULL) {

	for(auto it = ix_barcodes.begin(); it != ix_barcodes.end(); it++)
		this->_sample_map[it->first] = it->second.size() > 0 ? it->second[0] : "";
	post_init(reverse_complement);
}

unordered_set<int>* Barcode::get_set_sizes(){
  unordered_set<int>* set_sizes = new unordered_set<int>();
  auto it = allowed_lengths.find(length);
  if (it != allowed_lengths.end()){
      auto it_96 = lengths_96_barcodes.find(length);
      if(it_96 != lengths_96_barcodes.end())
          set_sizes->insert(96);
      auto it_384 = lengths_384_barcodes.find(length);
      if(it_384 != lengths_384_barcodes.end())
          set_sizes->insert(384);
  }
  return set_sizes;
}


void Barcode::post_init(bool reverse_complement){
    // if the barcode is reverse complement add _rc to the name
    // we need this later to load the correction maps
    if (reverse_complement)
        this->Name += "_rc";
    this->check_length();
}

void Barcode::check_length(){
// barcodes for each type need to be of the same length. This is as the
// sequences of 8 nt long barcodes are contained within 10 nt long barcodes,
// and the ones of 10 in 12. if we dont enforce the same length per barcode there
// might be a possibility we cant tell with certainty to which sample a
// barcoded read belongs.
	unordered_set<int> observed_lengths;
	for(auto it = this->_sample_map.begin(); it != this->_sample_map.end(); it++){
		string barcode = it->first;
		string sample_name = it->second;
		if (barcode.compare("") != 0){
			auto itl = this->allowed_lengths.find((int)barcode.length());
			if(itl != this->allowed_lengths.end()){
				observed_lengths.insert(barcode.length());
				if(observed_lengths.size() > 1){
					string message = string_format("%s barcodes with a different length "\
                                 "have been observed for %s. Barcodes"\
                                 " need to have the same length for all "\
                                 "samples.\nObserved barcode length:"\
                                 " %ld \nPrevious observed length:"\
                                 " %ld",sample_name.c_str(),sample_name.c_str(),barcode.length(), this->length);
					throw(std::runtime_error(message));
				}
				this->length = (int)barcode.length();
			}
			else{
				string tmplengths = "";
				for(auto ita = this->allowed_lengths.begin(); ita != this->allowed_lengths.end(); ita++)
					tmplengths += to_string(*ita) + ", ";
				string message = string_format("%s barcodes are %ld nt long.\n"\
                        "%s barcodes are only allowed to be"\
                        " %s nt long.", this->Name.c_str(), barcode.length(), this->Name.c_str(), tmplengths.substr(0,tmplengths.length()-1).c_str());
				throw(std::runtime_error(message));
			}

		}
	}
}

/**
 * Reads a tsv barcodes file, converts them to byes and returns them as a dict.
 These dicts are used to map erroneous barcodes to their corrected version.
 Args:
 barcode_type (str): A path pointing to a error correction map file.
 barcode_length (int): A path pointing to a error correction map file.

 Return:
 dict (dict): correction_map : <b'erroneous_barcode', b'corrected_barcode'>
 */
void Barcode::load_correction_map(string relative_exepath) {
	fprintf(stdout, "Trying to find the appropriate barcode set for %s...\n",
			this->Name.c_str());

	// When there are no barcodes specified, there is nothing to correct.
	if (this->length == 0) {
		fprintf(stdout, "No barcodes have been specified for %s.\n",
				this->Name.c_str());
		this->correction_map = new unordered_map<string, string>(); //.clear();
		this->correction_map->insert( { "", "" });
		return;
		//barcode.correction_map = NULL; //unordered_map<string,string> {"None": "None"};
		//return NULL;
	}
	unordered_set<int> *tmp_sizes = this->get_set_sizes();
	vector<int> sizes(tmp_sizes->begin(), tmp_sizes->end());
	delete tmp_sizes;
	std::sort(sizes.begin(), sizes.end());
	bool drop_none = true;
	int set_size;

	//size_t path_buffer = 10000;
	//char *buffer = (char*)malloc(sizeof(char)*path_buffer);
	//char * cwd = GETCWD(buffer, path_buffer);
	int index = relative_exepath.find_last_of(PATH_SEP);

	// try package path
	string package_str = relative_exepath.substr(0, index) + PATH_SEP + ".."
			+ PATH_SEP + "misc" + PATH_SEP + "barcodes" + PATH_SEP
			+ this->Name;
	if (!Parser::PathExists(package_str)) {
		// try install path
		package_str = relative_exepath.substr(0, index) + PATH_SEP + ".."
				+ PATH_SEP + "share" + PATH_SEP + "idemuxCPP" + PATH_SEP
				+ "barcodes" + PATH_SEP + this->Name;
		if (!Parser::PathExists(package_str)) {
			throw(runtime_error(
					string_format("Error: the path %s is not available!\n",
							package_str.c_str())));
		}
	}
	std::cout << package_str << std::endl;

	//free(cwd);

	for (auto it = sizes.begin(); it != sizes.end(); it++) {
		set_size = (*it);
		string file_str = "base_mapping_b" + to_string(set_size) + "_l"
				+ to_string(this->length) + ".tsv";

		unordered_map<string, string> *corr_map = Parser::get_map_from_resource(
				package_str, file_str);
		//_barcode_set = set(corr_map.values());
		unordered_set<string> *_barcodes_given = this->get_used_codes(
				drop_none);
		bool is_contained = true;
		for (auto it_test = _barcodes_given->begin();
				it_test != _barcodes_given->end(); it_test++) {
			auto it_contained = corr_map->find(*it_test);
			if (it_contained == corr_map->end()) {
				is_contained = false;
				break;
			}
		}
		delete _barcodes_given;
		if (is_contained) {
			printf(
					"Correct set found. Used set is %d barcodes with %d nt length.\n",
					set_size, this->length);
			//this->correction_map.clear();
			//for (auto it = corr_map->begin(); it != corr_map->end(); it++)
			//	this->correction_map[it->first] = it->second;
			this->correction_map = corr_map;
			//delete corr_map;
			return;
		}
		delete corr_map;
		/*
		 if _barcodes_given <= _barcode_set:
		 log.info(f"Correct set found. Used set is {set_size} barcodes with "
		 f"{this->length} nt length.")
		 barcode.correction_map = corr_map
		 return barcode;
		 */
	}
	if (this->length != 6) {
		printf(
				"Warning: No fitting Lexogen barcode set found for %s. No "
						"error correction will take place for this barcode. Are you using "
						"valid Lexogen barodes?\n", this->Name.c_str());
	}
	unordered_set<string> *_bc_list = this->used_codes();
	this->correction_map = new unordered_map<string, string>(); //.clear();
	for (auto it = _bc_list->begin(); it != _bc_list->end(); it++) {
		this->correction_map->insert( { *it, *it });
	}
	delete _bc_list;

	/*
	 if barcode.length != 6:
	 log.warning(f"No fitting Lexogen barcode set found for {barcode.name}. No "
	 f"error correction will take place for this barcode. Are you using "
	 f"valid Lexogen barodes?")

	 _bc_list = list(barcode.used_codes)
	 barcode.correction_map = dict(zip(_bc_list, _bc_list))
	 return barcode
	 */
}


Barcode::~Barcode() {
 	if(this->correction_map)
 		delete this->correction_map;
}

