#ifndef PAIREDREADER_H_
#define PAIREDREADER_H_

#include <string>
#include <vector>
#include "FastqReader.h"

using namespace std;

class PairedReader {
public:
	PairedReader(string fastqgz_1, string fastqgz_2);
	PairedReader(string bam_file_alternating);
	std::vector<std::pair<fq_read*, fq_read*>>* next_reads(size_t max_size);
	std::vector<std::pair<fq_read*, fq_read*>>* next_reads2(size_t max_size, int reading_threads = 2);
	virtual ~PairedReader();
private:
	IFastqReader *Reader1;
	IFastqReader *Reader2;
	bool IsBam;
};

#endif /* PAIREDREADER_H_ */
