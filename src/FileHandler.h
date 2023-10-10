#ifndef FILEHANDLER_H_
#define FILEHANDLER_H_
#include <unordered_map>
#include <string>
#include <cmath>
#include "ZipFastqWriter.h"

using namespace std;


class FileHandler {
public:
  FileHandler(unordered_map<string,string>& barcode_file_map, string output_folder, bool paired_output = false);
  void init_all_file_handles();
  void init_all_file_handles_single_end();
  std::pair<ZipFastqWriter*,ZipFastqWriter*>* get_file_handles(string &barcode);
  string get_sample_name(std::pair<ZipFastqWriter*,ZipFastqWriter*>* file_pair_key);
  virtual ~FileHandler();
private:
  unordered_map<string,string>* Barcode_file_map;
  string Output_folder;
  // map barcode to file handle
  unordered_map<string,std::pair<ZipFastqWriter*,ZipFastqWriter*>*> Fastq_handler;
  // map file handle to sample name
  unordered_map<std::pair<ZipFastqWriter*,ZipFastqWriter*>*, string> Map_filehandle_sample_name;
};

#endif /* FILEHANDLER_H_ */

