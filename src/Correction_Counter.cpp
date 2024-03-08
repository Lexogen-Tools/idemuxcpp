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
                                       unordered_set<string>  *i1codes)
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

    auto itcounts = Counts_combination.find(index);
    if (itcounts != Counts_combination.end()) {
      itcounts->second.all  += 1;
      itcounts->second.i1   += corrected_i1 ? 1 : 0;
      itcounts->second.i5   += corrected_i5 ? 1 : 0;
      itcounts->second.i7   += corrected_i7 ? 1 : 0;
    } else {
      Counts_combination[index].all = 1;
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
Correction_Counter::write_correction(string                         filename,
                                     unordered_map<string, string> &barcode_sample_map)
{
  ofstream        csvfile(filename, std::ofstream::out);

  string          csv_header =
    "sample_name,i7,i5,i1,combination_corrected,i7_corrected,i5_corrected,i1_corrected";

  csvfile << csv_header << std::endl;
  string          i7, i5, i1;
  barcode_counts  undetermined_counts = { 0, 0, 0, 0 };

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
    auto    it_combi  = barcode_sample_map.find(barcodes);
    if (it_combi != barcode_sample_map.end()) {
      string sample_name = it_combi->second;
      csvfile << sample_name << "," << i7 << "," << i5 << "," << i1 << "," << it->second.all <<
      "," << it->second.i7 << "," << it->second.i5 << "," << it->second.i1 << std::endl;
    } else {
      undetermined_counts.all += it->second.all;
      undetermined_counts.i7  += it->second.i7;
      undetermined_counts.i5  += it->second.i5;
      undetermined_counts.i1  += it->second.i1;
    }
  }
  if (undetermined_counts.all != 0) {
    string sample_name = "undetermined";
    csvfile << sample_name << ",,,," << undetermined_counts.all << "," << undetermined_counts.i7 <<
    "," << undetermined_counts.i5 << "," << undetermined_counts.i1 << std::endl;
  }

  csvfile.close();
  fprintf(stdout, "Barcode combination corrections summary saved to %s\n",
          filename.c_str());
}
