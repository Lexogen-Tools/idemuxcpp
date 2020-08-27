/*
 * ZipFastqWriter.cpp
 *
 *  Created on: Aug 26, 2020
 *      Author: gentzian
 */

#include <string>
#include <zlib.h>
#include "ZipFastqReader.h"
#include "ZipFastqWriter.h"


ZipFastqWriter::ZipFastqWriter(string outfile) : OutputFile(outfile) {
	// TODO Auto-generated constructor stub
	GZ_file_handle = gzopen(outfile.c_str(),"wb4");
	if(GZ_file_handle==NULL){ //throw std::runtime_error(strerror(errno));
		fprintf(stderr, "Error: could not open gz file for writing!\n");
		exit(EXIT_FAILURE);
	}
}

int ZipFastqWriter::write(const char* data, size_t len) {
    return gzwrite(GZ_file_handle, data, len);
}

void ZipFastqWriter::write_read(fq_read* r) {
		string read =  string(r->Seq_ID);
		read.append("\n");
		read.append(r->Sequence);
		read.append("\n");
		read.append(r->Plus_ID);
		read.append("\n");
		read.append(r->QualityCode);
		read.append("\n");
		auto res = this->write(read.data(), read.size());
		if(res==0) throw std::runtime_error("Error: could not write read from list to file!");
}

void ZipFastqWriter::write_read_list(fq_read** reads) {
	for(int i = 0; reads[i] != NULL; i++){
		fq_read *r = reads[i];
		string read =  string(r->Seq_ID);
		read.append("\n");
		read.append(r->Sequence);
		read.append("\n");
		read.append(r->Plus_ID);
		read.append("\n");
		read.append(r->QualityCode);
		read.append("\n");
		auto res = this->write(read.data(), read.size());
		if(res==0) throw std::runtime_error("Error: could not write read from list to file!");
	}
}

ZipFastqWriter::~ZipFastqWriter() {
	// TODO Auto-generated destructor stub
	gzclose(GZ_file_handle);
}

