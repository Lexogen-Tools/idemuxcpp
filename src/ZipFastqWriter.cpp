#include <string>
#include <zlib.h>
#include "ZipFastqReader.h"
#include "ZipFastqWriter.h"

ZipFastqWriter::ZipFastqWriter(string outfile) : OutputFile(outfile) {
	GZ_file_handle = gzopen(outfile.c_str(),"wb4");
	if(GZ_file_handle==NULL){
		string message = "Error: could not open gz file for writing!\n";
		throw(runtime_error(message));
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
		auto res = this->write(read.c_str(), read.length());
		if(res==0) throw std::runtime_error("Error: could not write read from list to file!");
}

void ZipFastqWriter::write_read_list(fq_read** reads) {
	for(int i = 0; reads[i] != NULL; i++){
		this->write_read(reads[i]);
	}
}

ZipFastqWriter::~ZipFastqWriter() {
	gzflush(GZ_file_handle,1);
	if (gzclose(GZ_file_handle) != Z_OK){
		std::cerr << "failed gzclose" << std::endl;
	}
}

