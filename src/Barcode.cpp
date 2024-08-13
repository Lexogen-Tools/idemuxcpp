#include <vector>
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <istream>
#include <fstream>
#include "config.h"
#include "helper.h"
#include "Parser.h"
#include "Barcode.h"
#include <boost/filesystem.hpp>
#include <set>
using namespace boost::filesystem;

using namespace std;

Barcode::Barcode(string                                       barcode_type,
                 unordered_map<string, std::vector<string> > &ix_barcodes,
                 bool                                         reverse_complement,
                 bool                                         auto_detect) :
  _Name(barcode_type), reverse_complement(reverse_complement), auto_detect(auto_detect),
  allowed_lengths({
  6, 8, 10, 12 }), lengths_96_barcodes({ 8, 10, 12 }), lengths_384_barcodes(
    { 10, 12 })
{
  for (auto it = ix_barcodes.begin(); it != ix_barcodes.end(); it++)
    this->_sample_map[it->first] =
      it->second.size() > 0 ? it->second[0] : "";

  unordered_set<int>* observed_lengths = this->read_used_lengths();
  this->Lengths.assign(observed_lengths->begin(), observed_lengths->end());
  delete observed_lengths;
}

string
Barcode::get_name(bool reverse)
{
  // if the barcode is reverse complement add _rc to the name
  // we need this later to load the correction maps
  if (this->reverse_complement != reverse)
    return this->_Name + "_rc";
  else
    return this->_Name;
}

unordered_set<int>*
Barcode::read_used_lengths()
{
  // barcodes for each type need to be of the same length. This is as the
  // sequences of 8 nt long barcodes are contained within 10 nt long barcodes,
  // and the ones of 10 in 12. if we dont enforce the same length per barcode there
  // might be a possibility we cant tell with certainty to which sample a
  // barcoded read belongs.
  unordered_set<int>* observed_lengths = new unordered_set<int>();

  for (auto it = this->_sample_map.begin(); it != this->_sample_map.end();
       it++) {
    string  barcode     = it->first;
    string  sample_name = it->second;
    if (barcode.compare("") != 0) {
      auto itl = this->allowed_lengths.find((int)barcode.length());
      if (itl == this->allowed_lengths.end()) {
        string  tmplengths = "";
        for (auto ita = this->allowed_lengths.begin();
             ita != this->allowed_lengths.end(); ita++)
          tmplengths += to_string(*ita) + ", ";
        string  message = string_format("%s barcodes are %ld nt long.\n"
                                        "%s barcodes are only allowed to be"
                                        " %s nt long. No error correction will be used.",
                                        this->get_name().c_str(),
                                        barcode.length(),
                                        this->get_name().c_str(),
                                        tmplengths.substr(0, tmplengths.length() - 1).c_str());
        //throw(std::runtime_error(message));
        printf("Warning: %s\n", message.c_str());
      }

      observed_lengths->insert(barcode.length());
    }
  }
  //this->Lengths.assign(observed_lengths.begin(), observed_lengths.end());
  return observed_lengths;
}


/**
 * Reads a tsv barcodes file, converts them to byes and returns them as a dict.
 * These dicts are used to map erroneous barcodes to their corrected version.
 * Args:
 * barcode_type (str): A path pointing to a error correction map file.
 * barcode_length (int): A path pointing to a error correction map file.
 *
 * Return:
 * dict (dict): correction_map : <b'erroneous_barcode', b'corrected_barcode'>
 */
void
Barcode::load_correction_map(string relative_exepath,
                             string correction_maps_path)
{
  fprintf(stdout, "Trying to find the appropriate barcode set for %s...\n",
          this->get_name().c_str());

  // When there are no barcodes specified, there is nothing to correct.
  if (this->Lengths.size() == 0) {
    fprintf(stdout, "No barcodes have been specified for %s.\n",
            this->get_name().c_str());
    this->Correction_map.clear();
    return;
  }

  bool                drop_none = true;
  string              barcode_path_prefix;

  if (!correction_maps_path.empty()) {
    barcode_path_prefix = correction_maps_path;
    if (!Parser::PathExists(barcode_path_prefix)) {
      throw(runtime_error(
              string_format("Error: the path %s is not available!\n",
                            barcode_path_prefix.c_str())));
    }
  } else {
    int index = relative_exepath.find_last_of(PATH_SEP);
    // try package path
    barcode_path_prefix = relative_exepath.substr(0, index) + PATH_SEP + ".."
                          + PATH_SEP + "misc" + PATH_SEP + "barcodes";
    if (!Parser::PathExists(barcode_path_prefix)) {
      // try install path
      barcode_path_prefix = relative_exepath.substr(0, index) + PATH_SEP + ".."
                            + PATH_SEP + "share" + PATH_SEP + PACKAGE_NAME + PATH_SEP
                            + "barcodes";
      if (!Parser::PathExists(barcode_path_prefix)) {
        // try /usr/share install path (/usr/local/share/idemuxcpp/barcodes/)
        barcode_path_prefix = string("") + PATH_SEP + "usr" + PATH_SEP + "local" + PATH_SEP
                              + PATH_SEP + "share" + PATH_SEP + PACKAGE_NAME + PATH_SEP
                              + "barcodes";
        if (!Parser::PathExists(barcode_path_prefix)) {
          // try /usr/share install path (/usr/share/idemuxcpp/barcodes/)
          barcode_path_prefix = PATH_SEP + string("usr")
                                + PATH_SEP + string("share") + PATH_SEP + PACKAGE_NAME + PATH_SEP
                                + string("barcodes");
          if (!Parser::PathExists(barcode_path_prefix)) {
            throw(runtime_error(
                    string_format("Error: the path %s is not available!\n",
                                  barcode_path_prefix.c_str())));
          }
        }
      }
    }
  }

  std::cout << barcode_path_prefix << std::endl;


  // sort codes according to length.
  unordered_set<string>                   *barcodes_given = this->get_used_codes(drop_none);
  string                                  code;
  unordered_map<size_t, vector<string> >  length_and_codes;

  for (auto it = barcodes_given->begin(); it != barcodes_given->end(); it++) {
    code = *it;
    length_and_codes[code.length()].push_back(code);
  }

  // load all barcodes for all set sizes and barcode lengths (that were defined in the sample sheet).
  this->Correction_map.clear();
  int   length;
  bool  contains_allowed_length = false;

  for (auto itl = this->Lengths.begin(); itl != this->Lengths.end(); itl++) {
    length = *itl;
    if (this->allowed_lengths.find(length) == this->allowed_lengths.end() || length == 6) {
      //do one to one mapping (there is no mapping table for this.)
      printf(
        "Barcodes with length %d found in dataset. These codes will be demultiplexed, but without error correction.\n",
        length);
      for (auto it_code_6 = length_and_codes[length].begin();
           it_code_6 != length_and_codes[length].end(); it_code_6++)
        this->Correction_map.insert({ *it_code_6, *it_code_6 });
      continue;
    } else {
      contains_allowed_length = true;
    }

    bool  reverse = false;

    /* Take smallest map, where all barcodes fit in:
     * 1. Read all maps for same length.
     * 2. Order by maps key set size.
     * 3. Take first, or error if second best set has same size.
    */
    unordered_map<set<string>*, unordered_map<string, string>*> bc_set_to_corr_map;

    //run twice if auto_detect is enabled to also check the rc
    for (int j = 0; j <= this->auto_detect; j++) {
      string dir_path = barcode_path_prefix + PATH_SEP + this->get_name(reverse) + PATH_SEP + "l" +
                        to_string(length);
      for (directory_entry& entry : directory_iterator(dir_path)) {
        unordered_map<string, string> *corr_map     = new unordered_map<string, string>();
        set<string>                   *barcode_set  = new set<string>();
        Parser::get_map_from_resource(entry.path().string(), corr_map, barcode_set);

        //test if all codes are contained in the map for a certain length.
        size_t                        n_codes_contained = 0;
        for (auto it_test = length_and_codes[length].begin();
             it_test != length_and_codes[length].end(); it_test++) {
          auto it_contained = barcode_set->find(*it_test);
          if (it_contained == barcode_set->end())
            break;

          n_codes_contained++;
        }

        if (n_codes_contained == length_and_codes[length].size()) {
            bc_set_to_corr_map[barcode_set] = corr_map;
            printf(
	    "Info: Matching set found (%ld matching codes). Used set is %ld barcodes with %d nt length.\n",
	    length_and_codes[length].size(),
	    barcode_set->size(),
	    length);
          //all codes for this length are contained in the set for the given size.
          // The real read set size is not always the set size from the file name. Thus, we don't break and collect all matching tables.
        }
        else{
            delete corr_map;
            delete barcode_set;
        }
      }
      //revert for the case that this loop is running twice
      reverse = !reverse;
    }


    vector<size_t> set_sizes;
    for(auto it = bc_set_to_corr_map.begin(); it != bc_set_to_corr_map.end(); it++){
        size_t set_size = it->first->size();
        set_sizes.push_back(set_size);
    }
    std::sort(set_sizes.begin(), set_sizes.end());

    if(set_sizes.size() > 0){
	size_t min_set_size = set_sizes[0];
	if(set_sizes.size() >= 2 && set_sizes[1] == min_set_size){
	    string message = "Error: found 2 barcode correction map of the same size (%s nt) that contain all given barcodes! Please use more barcodes or make the correction table unique!";
	    message = string_format(message, to_string(min_set_size));
	    throw(runtime_error(message));
	}

	set<string>* min_barcode_set = NULL;

	for(auto it = bc_set_to_corr_map.begin(); it != bc_set_to_corr_map.end(); it++){
	      size_t set_size = it->first->size();
	      if(set_size == min_set_size){
		  min_barcode_set = it->first;
		  break;
	      }
	}

	if(min_barcode_set != NULL){
	    length = min_barcode_set->begin()->length();
	    unordered_map<string, string> *corr_map = bc_set_to_corr_map[min_barcode_set];
	    printf(
	      "Correct set found (%ld matching codes). Used set is %ld barcodes with %d nt length.\n",
	      length_and_codes[length].size(),
	      min_barcode_set->size(),
	      length);
	    for (auto it_codes = corr_map->begin(); it_codes != corr_map->end(); it_codes++) {
	      auto it_found_code = this->Correction_map.find(it_codes->first);
	      if (it_found_code == this->Correction_map.end()) {
		this->Correction_map[it_codes->first] = it_codes->second;
	      }
	    }
	}
	else{
	    string message = "Error: could not find a valid correction set!";
	    throw(runtime_error(message));
	}
    }

    for(auto it = bc_set_to_corr_map.begin(); it != bc_set_to_corr_map.end(); it++){
            delete it->first;
            delete it->second;
    }
  }

  if (!contains_allowed_length) {
    printf(
      "Warning: No fitting Lexogen barcode set found for %s. No "
      "error correction will take place for this barcode. Are you using "
      "valid Lexogen barcodes?\n", this->get_name().c_str());
  }

  if (this->Correction_map.size() == 0) {
    for (auto it = barcodes_given->begin(); it != barcodes_given->end(); it++)
      this->Correction_map.insert({ *it, *it });
  }

  delete barcodes_given;
}


void
Barcode::create_one_to_one_map()
{
  if (this->Lengths.size() > 0) {
    unordered_set<string> *_bc_list = this->used_codes();
    this->Correction_map.clear();
    for (auto it = _bc_list->begin(); it != _bc_list->end(); it++)
      this->Correction_map.insert({ *it, *it });
    delete _bc_list;
  } else {
    // When there are no barcodes specified, there is nothing to correct.
    fprintf(stdout, "No barcodes have been specified for %s.\n",
            this->get_name().c_str());
    this->Correction_map.clear();
    this->Correction_map.insert({ "", "" });
  }
}


Barcode::~Barcode()
{
  //if (this->Correction_map)
  //	delete this->Correction_map;
}
