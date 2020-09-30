/*
 * FastqReader.cpp
 *
 *  Created on: Oct 29, 2020
 *      Author: gentzian
 */

#include "FastqReader.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

using namespace std;

FastqReader::FastqReader(string filename) {
	FastqFileHandle.open(filename, std::ifstream::in);
}

fq_read* FastqReader::next_read() {
	std::string line;
	fq_read *read = new fq_read();
	short i = 0;
	while (!FastqFileHandle.eof() && i < 4) {
		std::getline(FastqFileHandle, line);
		if (FastqFileHandle.bad() || FastqFileHandle.fail()) {
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

void FastqReader::close() {
	FastqFileHandle.close();
}

FastqReader::~FastqReader() {
	this->close();
}

