#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <zlib.h>
#include <stdexcept>
#include "FastqReader.h"
#include "ZipFastqWriter.h"

/**
 * Compress memory via zlib and include header (equivalent to block-gzip).
 * Inspired by https://stackoverflow.com/questions/4538586/how-to-compress-a-buffer-with-zlib and https://www.zlib.net/zlib_how.html
 *
 * @param in_data - arbitrary array of input data.
 * @param in_data_size - length of input data.
 * @param out_data - already allocated (new()-ed) output data.
 * @param compression_level - zlib compression level (value between 0 and 9).
 */
int ZipFastqWriter::compress_memory(uint8_t *in_data, size_t in_data_size, std::vector<uint8_t> *out_data, int compression_level)
{
	int result = 0;
	std::vector<uint8_t> *buffer = out_data;

	const size_t BUFSIZE = 128 * 1024;
	uint8_t temp_buffer[BUFSIZE];

	z_stream strm;
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	strm.next_in = reinterpret_cast<uint8_t *>(in_data);
	strm.avail_in = in_data_size;
	strm.next_out = temp_buffer;
	strm.avail_out = BUFSIZE;

	int windowBits = 15; // default value according to zlib manual.
	windowBits += 16; // 16 means that it should write the header and crc32.

	if(Z_OK != deflateInit2(&strm, compression_level,Z_DEFLATED, windowBits, 9, Z_DEFAULT_STRATEGY)){
		std::cerr << "Error: no gzip init!" << std::endl;
	}

	gz_header header {0};
	header.name = Z_NULL;
	header.comment = Z_NULL;
	header.extra = Z_NULL;
	if(Z_OK != deflateSetHeader(&strm, &header)){
		std::cerr << "Error: no zlib deflateSetHeader!" << std::endl;
	}
	while (strm.avail_in != 0)
	{
		int res = deflate(&strm, Z_NO_FLUSH);
		if(res == Z_OK){
			//std::cerr << "deflate Z_OK" << std::endl;
		}
		if (strm.avail_out == 0)
		{
			buffer->insert(buffer->end(), temp_buffer, temp_buffer + BUFSIZE);
			strm.next_out = temp_buffer;
			strm.avail_out = BUFSIZE;
		}
	}

	int deflate_res = Z_OK;
	while (deflate_res == Z_OK)
	{
		if (strm.avail_out == 0)
		{
			buffer->insert(buffer->end(), temp_buffer, temp_buffer + BUFSIZE);
			strm.next_out = temp_buffer;
			strm.avail_out = BUFSIZE;
		}
		deflate_res = deflate(&strm, Z_FINISH);
	}

	if(deflate_res == Z_STREAM_END){
		//std::cerr << "gz stream end" << std::endl;
	}
	buffer->insert(buffer->end(), temp_buffer, temp_buffer + BUFSIZE - strm.avail_out);
	deflateEnd(&strm);

	return result;
}


ZipFastqWriter::ZipFastqWriter(string outfile) : OutputFile(outfile) {
	File_handle = fopen(outfile.c_str(),"w");
	if(File_handle==NULL){
		string message = "Error: could not open gz file for writing!\n";
		throw(runtime_error(message));
	}

	fflush(File_handle);
	if (fclose(File_handle) != 0){
		std::cerr << "failed gzclose" << std::endl;
	}
}


int ZipFastqWriter::write(std::vector<uint8_t> *out_data){
	int result = 0;
	File_handle = fopen(OutputFile.c_str(), "a");
	if(File_handle==NULL){
		string message = "Error: could not open gz file for writing!\n";
		throw(runtime_error(message));
	}
	size_t compressed_size = out_data->size();
	uint8_t* out_array = (uint8_t*)malloc(sizeof(uint8_t)*compressed_size);
	copy(out_data->begin(), out_data->end(), out_array);
	fwrite((void*)out_array, sizeof(uint8_t), compressed_size, File_handle);
	free(out_array);
	if (fclose(File_handle) != 0){
		std::cerr << "failed gzclose" << std::endl;
		result = 1;
	}
	return result;
}

ZipFastqWriter::~ZipFastqWriter() {
}

