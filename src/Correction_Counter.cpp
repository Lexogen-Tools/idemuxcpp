#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <stdexcept>
#include "Correction_Counter.h"

using namespace std;

Correction_Counter::Correction_Counter(unordered_set<string>  *i7codes,
                                       unordered_set<string>  *i5codes,
                                       unordered_set<string>  *i1codes,
				       const unordered_map<string, string>   *barcode_sample_map) : Barcode_sample_map(barcode_sample_map)
{
  uint16_t i = 0;

  for (auto it = i7codes->begin(); it != i7codes->end(); it++) {
    if ((*it).compare("") != 0) {
      I7_codes.push_back(*it);
      I7_indices[*it] = i;
      i++;
    }
  }
  i = 0;
  for (auto it = i5codes->begin(); it != i5codes->end(); it++) {
    if ((*it).compare("") != 0) {
      I5_codes.push_back(*it);
      I5_indices[*it] = i;
      i++;
    }
  }
  i = 0;
  for (auto it = i1codes->begin(); it != i1codes->end(); it++) {
    if ((*it).compare("") != 0) {
      I1_codes.push_back(*it);
      I1_indices[*it] = i;
      i++;
    }
  }
};

void
Correction_Counter::count_correction(string i7,
                                     string i5,
                                     string i1,
                                     bool   corrected_i7,
                                     bool   corrected_i5,
                                     bool   corrected_i1)
{
  if (corrected_i7 || corrected_i5 || corrected_i1) {
    auto      iti7  = I7_indices.find(i7);
    auto      iti5  = I5_indices.find(i5);
    auto      iti1  = I1_indices.find(i1);
    uint64_t  index = 0;
    if (iti7 != I7_indices.end())
      index += iti7->second + 1;

    index = index << 16;
    if (iti5 != I5_indices.end())
      index += iti5->second + 1;

    index = index << 16;
    if (iti1 != I1_indices.end())
      index += iti1->second + 1;

    corrected_i7 = corrected_i7 && i7 != "" && iti7 != I7_indices.end();
    corrected_i5 = corrected_i5 && i5 != "" && iti5 != I5_indices.end();
    corrected_i1 = corrected_i1 && i1 != "" && iti1 != I1_indices.end();
    string  barcodes  = i7 + '\n' + i5 + '\n' + i1;
    auto it_combi  = Barcode_sample_map->find(barcodes);
    bool is_wanted_combination = (it_combi != Barcode_sample_map->end());
    bool combination_corrected = (corrected_i7 || corrected_i5 || corrected_i1) && is_wanted_combination;

    auto itcounts = Counts_combination.find(index);
    if (itcounts != Counts_combination.end()) {
      itcounts->second.all  += combination_corrected ? 1 : 0;
      itcounts->second.i1   += corrected_i1 ? 1 : 0;
      itcounts->second.i5   += corrected_i5 ? 1 : 0;
      itcounts->second.i7   += corrected_i7 ? 1 : 0;
    } else {
      Counts_combination[index].all = combination_corrected ? 1 : 0;
      Counts_combination[index].i1  = corrected_i1 ? 1 : 0;
      Counts_combination[index].i5  = corrected_i5 ? 1 : 0;
      Counts_combination[index].i7  = corrected_i7 ? 1 : 0;
    }
  }
}


Correction_Counter&
Correction_Counter::operator+=(const Correction_Counter& toAdd)
{
  //check if barcode indices are the same.
  bool    b_same = true;

  size_t  i;

  for (i = 0; i < toAdd.I7_codes.size(); i++) {
    if (this->I7_codes[i] != toAdd.I7_codes[i]) {
      b_same = false;
      break;
    }
  }
  for (i = 0; i < toAdd.I5_codes.size(); i++) {
    if (this->I5_codes[i] != toAdd.I5_codes[i]) {
      b_same = false;
      break;
    }
  }
  for (i = 0; i < toAdd.I1_codes.size(); i++) {
    if (this->I1_codes[i] != toAdd.I1_codes[i]) {
      b_same = false;
      break;
    }
  }
  if (!b_same) {
    fprintf(stderr, "Error: cannot add correction counters with different barcode sets!\n");
    exit(EXIT_FAILURE);
  }

  for (auto it = toAdd.I7_indices.begin(); it != toAdd.I7_indices.end(); it++) {
    if (this->I7_indices.find(it->first) != this->I7_indices.end())
      this->I7_indices[it->first] += it->second;
    else
      this->I7_indices[it->first] = it->second;
  }

  for (auto it = toAdd.I5_indices.begin(); it != toAdd.I5_indices.end(); it++) {
    if (this->I5_indices.find(it->first) != this->I5_indices.end())
      this->I5_indices[it->first] += it->second;
    else
      this->I5_indices[it->first] = it->second;
  }

  for (auto it = toAdd.I1_indices.begin(); it != toAdd.I1_indices.end(); it++) {
    if (this->I1_indices.find(it->first) != this->I1_indices.end())
      this->I1_indices[it->first] += it->second;
    else
      this->I1_indices[it->first] = it->second;
  }

  for (auto it = toAdd.Counts_combination.begin(); it != toAdd.Counts_combination.end(); it++) {
    if (this->Counts_combination.find(it->first) != this->Counts_combination.end())
      this->Counts_combination[it->first] += it->second;
    else
      this->Counts_combination[it->first] = it->second;
  }
  return *this;
}


/*
 * each line contains the corrected code, if it was corrected.
 */
void
Correction_Counter::write_correction(string filename)
{
  ofstream        csvfile(filename, std::ofstream::out);

  string          csv_header =
    "sample_name,i7,i5,i1,combination_corrected,i7_corrected,i5_corrected,i1_corrected";

  csvfile << csv_header << std::endl;
  string          i7, i5, i1;
  barcode_counts  undetermined_counts = { 0, 0, 0, 0 };

  vector<std::pair<string, string>> sorted_sample_names;
  for(auto it = Barcode_sample_map->begin(); it != Barcode_sample_map->end(); it++){
      sorted_sample_names.push_back(std::pair<string,string>(it->second, it->first));
  }
  std::sort(sorted_sample_names.begin(), sorted_sample_names.end());

  unordered_map<string,barcode_counts> sample_name_counts;
  for (auto it = Counts_combination.begin(); it != Counts_combination.end(); it++) {
    uint64_t  code_idc  = it->first;
    uint16_t  ui_i1     = uint16_t(code_idc);
    uint16_t  ui_i5     = uint16_t(code_idc >> 16);
    uint16_t  ui_i7     = uint16_t(code_idc >> 32);
    if (ui_i1 > 0)
      i1 = I1_codes[ui_i1 - 1];
    else
      i1 = "";

    if (ui_i5 > 0)
      i5 = I5_codes[ui_i5 - 1];
    else
      i5 = "";

    if (ui_i7 > 0)
      i7 = I7_codes[ui_i7 - 1];
    else
      i7 = "";

    // if combi in samples, else add to undetermined.
    string  barcodes  = i7 + '\n' + i5 + '\n' + i1;
    auto    it_combi  = Barcode_sample_map->find(barcodes);
    if (it_combi != Barcode_sample_map->end()) {
      string sample_name = it_combi->second;
      sample_name_counts[sample_name] = it->second;
    }
  }
  for(auto it = sorted_sample_names.begin(); it != sorted_sample_names.end(); it++){
      string sample_name = it->first;
      if(sample_name.compare("undetermined") == 0)
      	  continue;
      string str_codes = it->second;
      int first_index = str_codes.find('\n');
      if(first_index < 0){
	  i7 = i5 = i1 = "";
      }
      else{
	  i7 = str_codes.substr(0, first_index);
	  string i5_i1 = str_codes.substr(first_index+1);
	  int second_index = i5_i1.find('\n');
	  if(second_index < 0){
	      i5 = str_codes;
	      i1 = "";
	  }
	  else{
	      i5 = i5_i1.substr(0, second_index);
	      i1 = i5_i1.substr(second_index+1);
	  }
      }
      auto it_bc = sample_name_counts.find(sample_name);
      barcode_counts counts;
      if(it_bc != sample_name_counts.end()){
	  counts = it_bc->second;
      }
      else{
	  counts = { 0, 0, 0, 0 };
      }
      csvfile << sample_name << "," << i7 << "," << i5 << "," << i1 << "," << counts.all <<
	      "," << counts.i7 << "," << counts.i5 << "," << counts.i1 << std::endl;
  }

  csvfile.close();
  fprintf(stdout, "Barcode combination corrections summary saved to %s\n",
          filename.c_str());
}
