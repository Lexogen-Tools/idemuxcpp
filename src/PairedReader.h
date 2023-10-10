#ifndef PAIREDREADER_H_
#define PAIREDREADER_H_

#include <string>
#include <vector>
#include "FastqReader.h"

using namespace std;

class PairedReader {
public:
	PairedReader(string fastqgz_1, string fastqgz_2, double queue_buffer_gb);
	PairedReader(string single_read_file, double queue_buffer_gb, bool is_interleaved);
	std::vector<std::pair<fq_read*, fq_read*>>* next_reads(size_t max_size);
	std::vector<std::pair<fq_read*, fq_read*>>* next_reads2(size_t max_size, int reading_threads = 2);
	virtual ~PairedReader();
private:
	void fill_reads_interleaved(IFastqReader* reader, size_t max_size, std::vector<std::pair<fq_read*, fq_read*>> *pairs);
	size_t Queue_buffer_bytes;
	bool IsInterleaved;
	IFastqReader *Reader1;
	IFastqReader *Reader2;
};

#endif /* PAIREDREADER_H_ */
