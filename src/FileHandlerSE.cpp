
#include <unordered_map>
#include <string>
#include <cmath>
#include "helper.h"
#include "FileHandlerSE.h"

FileHandlerSE::FileHandlerSE(unordered_map<string,string>& barcode_file_map, string output_folder, double size_writer_buffer_gb) {
	Barcode_file_map = &barcode_file_map;
	(*Barcode_file_map)["undetermined"] = "undetermined";
	Output_folder = output_folder;
	if (!utils::folder_exists(output_folder))
		utils::mkdir(output_folder.c_str());
	init_all_file_handles(size_writer_buffer_gb);
}

// return and open file handles
void FileHandlerSE::init_all_file_handles(double total_buffer_gb){
	string barcode, sample_name;
	size_t given_buffer_size_per_file_bytes = int(floor(((total_buffer_gb/(Barcode_file_map->size()*2))*pow(1024,3))));
	size_t buffer_size_bytes = max(size_t(2), given_buffer_size_per_file_bytes);

	for(auto it = Barcode_file_map->begin(); it != Barcode_file_map->end(); it++){
		barcode = it->first;
		sample_name = it->second;
		string mate1_path = Output_folder +"/"+ sample_name + ".fastq.gz";
		ZipFastqWriter *w1 = new ZipFastqWriter(mate1_path, buffer_size_bytes);
		Fastq_handler[barcode] = w1;
		Map_filehandle_sample_name[w1] = sample_name;
	}
}

ZipFastqWriter* FileHandlerSE::get_file_handles(string &barcode){
	auto it = Fastq_handler.find(barcode);
	if (it != Fastq_handler.end())
		return it->second;
	else
		return Fastq_handler["undetermined"];
	return NULL;
}

string FileHandlerSE::get_sample_name(ZipFastqWriter* file_key){
	auto it = Map_filehandle_sample_name.find(file_key);
	if (it != Map_filehandle_sample_name.end())
		return it->second;
	return "";
}

FileHandlerSE::~FileHandlerSE() {
	// clear file handles
	for(auto it = Fastq_handler.begin(); it != Fastq_handler.end(); it++){
		delete it->second;
	}
}

