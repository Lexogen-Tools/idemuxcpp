
#ifndef SRC_BAMREADER_H_
#define SRC_BAMREADER_H_

#include <fstream>
#include <api/BamReader.h>
#include "FastqReader.h"


using namespace std;

class BamReader : public IFastqReader{
public:
	BamReader(string filename);
	fq_read* next_read();
        std::pair<fq_read*, fq_read*> next_read_pair();
	void close();
	virtual ~BamReader();
private:
        fq_read* get_next_read_from_bam(string &out_i7, string &out_i5, int &out_read_type);
	BamTools::BamReader* FileHandle;
};

#endif /* SRC_BOOSTZIPREADER_H_ */
