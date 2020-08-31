#ifndef PAIREDREADER_H_
#define PAIREDREADER_H_

#include <string>
#include <vector>
#include "ZipFastqReader.h"

using namespace std;

class PairedReader {
public:
	PairedReader(string fastqgz_1, string fastqgz_2);
	std::vector<std::pair<fq_read*, fq_read*>>* next_reads(size_t max_size);
	std::vector<std::pair<fq_read*, fq_read*>>* next_reads2(size_t max_size, int reading_threads = 2);
	virtual ~PairedReader();
private:
	ZipFastqReader Reader1;
	ZipFastqReader Reader2;
};

#endif /* PAIREDREADER_H_ */
