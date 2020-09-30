
#ifndef SRC_BOOSTZIPREADER_H_
#define SRC_BOOSTZIPREADER_H_

#include <fstream>
#include <boost/iostreams/filtering_stream.hpp>
#include "FastqReader.h"


using namespace std;

class BoostZipReader : public IFastqReader{
public:
	BoostZipReader(string filename);
	fq_read* next_read();
	void close();
	virtual ~BoostZipReader();
private:
	fstream FileHandle;
	boost::iostreams::filtering_istream GzipFileHandle;

};

#endif /* SRC_BOOSTZIPREADER_H_ */
