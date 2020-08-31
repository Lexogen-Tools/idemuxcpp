#include "PairedReader.h"

PairedReader::PairedReader(string fastqgz_1, string fastqgz_2) :
		Reader1(fastqgz_1), Reader2(fastqgz_2) {
	// TODO Auto-generated constructor stub

}

std::vector<std::pair<fq_read*, fq_read*>>* PairedReader::next_reads(
		size_t max_size) {
	std::vector<std::pair<fq_read*, fq_read*>> *pairs = new std::vector<std::pair<fq_read*, fq_read*>>();
	for (size_t i = 0; i < max_size; i++) {
		fq_read *read1 = Reader1.next_read(NULL);
		fq_read *read2 = Reader2.next_read(NULL);
		if (read1 != NULL && read2 != NULL) {
			pairs->push_back(std::pair<fq_read*, fq_read*>(read1, read2));
		}
		else{
			delete read1;
			delete read2;
		}

	}
	return pairs;
}

PairedReader::~PairedReader() {
	// TODO Auto-generated destructor stub
	//this->Reader1.close();
	//this->Reader2.close();
}

