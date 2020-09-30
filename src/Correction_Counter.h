/*
 * Correction_Counter.h
 *
 *  Created on: Oct 5, 2020
 *      Author: gentzian
 */

#ifndef SRC_CORRECTION_COUNTER_H_
#define SRC_CORRECTION_COUNTER_H_

#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <stdexcept>

using namespace std;

enum barcode_type{
	i1,
	i5,
	i7
};

struct barcode_counts{
	// corrections within one combination.
	uint64_t i1;
	uint64_t i5;
	uint64_t i7;
	//counted the whole combination.
	uint64_t all;
};

class Correction_Counter{
public:
	Correction_Counter(unordered_set<string>* i7codes, unordered_set<string>* i5codes, unordered_set<string>* i1codes);

	/**
	 * Counts the corrections for the combination of the given indices.
	 * if a string in this combination has not been used, insert the empty string "".
	 */
	void count_correction(string i7, string i5, string i1, bool corrected_i7, bool corrected_i5, bool corrected_i1);

	/*
	 * each line contains the corrected code, if it was corrected.
	 */
	void write_correction(string filename, unordered_map<string, string> &barcode_sample_map);

private:
	vector<string> I7_codes;
	vector<string> I5_codes;
	vector<string> I1_codes;
	unordered_map<string, uint16_t> I7_indices;
	unordered_map<string, uint16_t> I5_indices;
	unordered_map<string, uint16_t> I1_indices;
	unordered_map<uint64_t, barcode_counts> Counts_combination;
};

#endif /* SRC_CORRECTION_COUNTER_H_ */
