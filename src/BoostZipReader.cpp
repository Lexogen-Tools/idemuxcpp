
#include <iostream>
#include <fstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include "FastqReader.h"
#include "BoostZipReader.h"

using namespace std;

BoostZipReader::BoostZipReader(string filename) {
	this->FileHandle.open(filename, std::ios_base::in | std::ios_base::binary);
	this->GzipFileHandle.push(boost::iostreams::gzip_decompressor());
	this->GzipFileHandle.push(this->FileHandle);
}

fq_read* BoostZipReader::next_read() {
	std::string line;
	fq_read *read = new fq_read();
	short i = 0;
	while (!this->GzipFileHandle.eof() && i < 4) {
		std::getline(this->GzipFileHandle, line);
		if (this->GzipFileHandle.bad() || this->GzipFileHandle.fail()) {
			//(runtime_error("Error: could not read from fastq file!"));
			break;
		}
		switch(i){
		case 0:
			read->Seq_ID = line;
			break;
		case 1:
			read->Sequence = line;
			break;
		case 2:
			read->Plus_ID = line;
			break;
		case 3:
			read->QualityCode = line;
			break;
		default:
			break;
		}
		i++;
	}
	if (i < 4){
		delete read;
		return NULL;
	}
	return read;
}

void BoostZipReader::close() {
	this->FileHandle.close();
}

BoostZipReader::~BoostZipReader() {
	this->close();
}

