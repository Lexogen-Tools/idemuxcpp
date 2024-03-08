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
  post_init();
}


unordered_set<int> *
Barcode::get_set_sizes()
{
  unordered_set<int>  *set_sizes = new unordered_set<int>();
  int                 length;

  for (size_t i = 0; i < this->Lengths.size(); i++) {
    length = this->Lengths[i];
    auto it = allowed_lengths.find(length);
    if (it != allowed_lengths.end()) {
      auto  it_96 = lengths_96_barcodes.find(length);
      if (it_96 != lengths_96_barcodes.end())
        set_sizes->insert(96);

      auto  it_384 = lengths_384_barcodes.find(length);
      if (it_384 != lengths_384_barcodes.end())
        set_sizes->insert(384);
    }
  }
  return set_sizes;
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


void
Barcode::post_init()
{
  this->check_length();
}


void
Barcode::check_length()
{
  // barcodes for each type need to be of the same length. This is as the
  // sequences of 8 nt long barcodes are contained within 10 nt long barcodes,
  // and the ones of 10 in 12. if we dont enforce the same length per barcode there
  // might be a possibility we cant tell with certainty to which sample a
  // barcoded read belongs.
  unordered_set<int> observed_lengths;

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

      observed_lengths.insert(barcode.length());
    }
  }
  this->Lengths.assign(observed_lengths.begin(), observed_lengths.end());
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
    //this->Correction_map.insert( { "", "" });
    return;
    //barcode.correction_map = NULL; //unordered_map<string,string> {"None": "None"};
    //return NULL;
  }

  unordered_set<int>  *tmp_sizes = this->get_set_sizes();
  vector<int>         sizes(tmp_sizes->begin(), tmp_sizes->end());

  delete tmp_sizes;
  std::sort(sizes.begin(), sizes.end());
  bool                drop_none = true;

  //size_t path_buffer = 10000;
  //char *buffer = (char*)malloc(sizeof(char)*path_buffer);
  //char * cwd = GETCWD(buffer, path_buffer);
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
  int   n_loaded_maps = 0;
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
    bool  found   = false;
    //run twice if auto_detect is enabled to also check the rc
    for (int j = 0; j <= this->auto_detect && found == false; j++) {
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
          printf(
            "Correct set found (%ld matching codes). Used set is %ld barcodes with %d nt length.\n",
            length_and_codes[length].size(),
            barcode_set->size(),
            length);
          for (auto it_codes = corr_map->begin(); it_codes != corr_map->end(); it_codes++) {
            auto it_found_code = this->Correction_map.find(it_codes->first);
            if (it_found_code == this->Correction_map.end()) {
              this->Correction_map[it_codes->first] = it_codes->second;
            } else if (it_found_code->second != it_codes->second) {
              delete corr_map;
              delete barcode_set;
              string message =
                "Error: the sample sheet contains barcodes, which cannot be uniquely mapped to your correction tables!\n"
                "The barcode %s is mapped to %s in the current correction table and to %s in the previous correction table!";
              message = string_format(message,
                                      it_codes->first.c_str(),
                                      it_codes->second.c_str(),
                                      it_found_code->second.c_str());
              throw(runtime_error(message));
            }
          }
          n_loaded_maps++;
          //all codes for this length are contained in the set for the given size --> break.
          found = true;
        }

        delete corr_map;
        delete barcode_set;
      }
      //revert for the case that this loop is running twice
      reverse = !reverse;
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
