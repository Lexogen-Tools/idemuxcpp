#include <string>
#include <iostream>
#include <zlib.h>
#include <stdexcept>
#include "FastqReader.h"
#include "ZipFastqWriter.h"

ZipFastqWriter::ZipFastqWriter(string outfile, size_t buffer_size_bytes) : OutputFile(outfile), Buffer("") {
	size_t max_size = size_t(-1);
	this->BufferSize = min(max(buffer_size_bytes,1ul),max_size-1);
	// increase string capacity
	this->Buffer.reserve(BufferSize);

	GZ_file_handle = gzopen(outfile.c_str(),"wb4");
	if(GZ_file_handle==NULL){
		string message = "Error: could not open gz file for writing!\n";
		throw(runtime_error(message));
	}

	gzflush(GZ_file_handle,1);
	if (gzclose(GZ_file_handle) != Z_OK){
		std::cerr << "failed gzclose" << std::endl;
	}
}

int ZipFastqWriter::write_fragments(const char* data, size_t len, size_t buffer_size){
	int result = 0;
	size_t max_gz_write = unsigned(-1) -1;
	size_t new_buffer_size = min(buffer_size, max_gz_write);

	GZ_file_handle = gzopen(OutputFile.c_str(),"ab4");
	if(GZ_file_handle==NULL){
		string message = "Error: could not open gz file for writing!\n";
		throw(runtime_error(message));
	}

	// data does not even fit into buffer --> directly write it.
	if(len > new_buffer_size){
		size_t offset;
		unsigned current_length;
		for(size_t i = 0; i < len; i+=new_buffer_size){
			offset = i;
			current_length = unsigned(min(new_buffer_size, len - offset));
			result = result | gzwrite(GZ_file_handle, &data[offset], current_length);
		}
	}
	else{
		gzwrite(GZ_file_handle, data, len);
	}

	gzflush(GZ_file_handle,1);
	if (gzclose(GZ_file_handle) != Z_OK){
		std::cerr << "failed gzclose" << std::endl;
	}
	return result;
}

int ZipFastqWriter::write(const char* data, size_t len) {
	int result = 0;
	size_t buffer_length = this->Buffer.length();
	// 1. only append to buffer and write in ZipFastqWriter::flush()
	if(buffer_length/2 + len/2 < BufferSize/2){
		this->Buffer.append(data);
		std::cout << "append" << std::endl;
	}
	else{
		// 2. Buffer would overflow --> write buffer, clear it and append.
		this->flush();
		buffer_length = this->Buffer.length();
		if(buffer_length/2 + len/2 < BufferSize/2){
			this->Buffer.append(data);
		}
		else{
			// 3. if buffer is to small for input data, write it directly.
			//if(len >= BufferSize){
			this->write_fragments(data, len, len);
		}
	}
	return result;
}

void ZipFastqWriter::write_read(fq_read* r) {
		string read = r->to_string();
		auto res = this->write(read.c_str(), read.length());
		if(res==0) throw std::runtime_error("Error: could not write read from list to file!");
}

void ZipFastqWriter::write_read_list(fq_read** reads) {
	for(int i = 0; reads[i] != NULL; i++){
		this->write_read(reads[i]);
	}
}

void ZipFastqWriter::flush(){
	size_t buffer_length = this->Buffer.length();
	if(buffer_length > 0){
		this->write_fragments(Buffer.c_str(), buffer_length, buffer_length);
		this->Buffer.clear();
	}
}

ZipFastqWriter::~ZipFastqWriter() {
	this->flush();
	/*
	gzflush(GZ_file_handle,1);
	if (gzclose(GZ_file_handle) != Z_OK){
		std::cerr << "failed gzclose" << std::endl;
	}
	*/
}

