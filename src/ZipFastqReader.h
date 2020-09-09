#ifndef ZIPFASTQREADER_H_
#define ZIPFASTQREADER_H_

#include <string>
#include <zlib.h>
#include <iostream>


using namespace std;

struct fq_read{
	string Seq_ID;
	string Sequence;
	string Plus_ID;
	string QualityCode;
	fq_read() : Seq_ID(""),Sequence(""),Plus_ID(""),QualityCode("") {};
	~fq_read(){

	};
	inline string to_string(){
		string read =  string(this->Seq_ID);
				read.append("\n");
				read.append(this->Sequence);
				read.append("\n");
				read.append(this->Plus_ID);
				read.append("\n");
				read.append(this->QualityCode);
				read.append("\n");
				return read;
	};
};

class ZipFastqReader {
public:
	ZipFastqReader(string filename);
	fq_read* next_read();
	void close();
	virtual ~ZipFastqReader();
private:
	gzFile GZ_filehandle;
	size_t BUFLEN;
	char* buf;
	char* offset;
	char* cur;
	char* end;
	int read_line_index;
};

#endif /* ZIPFASTQREADER_H_ */
