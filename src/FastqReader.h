/*
 * FastqReader.h
 *
 *  Created on: Oct 29, 2020
 *      Author: gentzian
 */

#ifndef SRC_FASTQREADER_H_
#define SRC_FASTQREADER_H_

#include <string>
#include <fstream>

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

class IFastqReader{
public:
	virtual fq_read* next_read() = 0;
	virtual std::pair<fq_read*, fq_read*> next_read_pair(){ return std::pair<fq_read*, fq_read*>({NULL, NULL}); };
	virtual ~IFastqReader(){};
};

class FastqReader : public IFastqReader{
public:
	FastqReader(string filename);
	fq_read* next_read();
	void close();
	virtual ~FastqReader();
private:
	std::ifstream FastqFileHandle;
};

#endif /* SRC_FASTQREADER_H_ */
