#ifndef SRC_FILEHANDLERSE_H_
#define SRC_FILEHANDLERSE_H_

#include <unordered_map>
#include <string>
#include <cmath>
#include "ZipFastqWriter.h"

using namespace std;

class FileHandlerSE {
public:
	FileHandlerSE(unordered_map<string, string> &barcode_file_map,
			string output_folder, double size_writer_buffer_gb);
	void init_all_file_handles(double total_buffer_gb);
	ZipFastqWriter* get_file_handles(string &barcode);
	string get_sample_name(ZipFastqWriter *file_pair_key);
	virtual ~FileHandlerSE();
private:
	unordered_map<string, string> *Barcode_file_map;
	string Output_folder;
	// map barcode to file handle
	unordered_map<string, ZipFastqWriter*> Fastq_handler;
	// map file handle to sample name
	unordered_map<ZipFastqWriter*, string> Map_filehandle_sample_name;

};

#endif /* SRC_FILEHANDLERSE_H_ */
