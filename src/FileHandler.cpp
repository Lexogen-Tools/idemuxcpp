#include "FileHandler.h"
#include <unordered_map>
#include <string>
#include <cmath>
#include <algorithm>
#include "helper.h"

FileHandler::FileHandler(unordered_map<string, string>& barcode_file_map,
                         string                         output_folder,
                         bool                           paired_output)
{
  Barcode_file_map                    = &barcode_file_map;
  (*Barcode_file_map)["undetermined"] = "undetermined";
  Output_folder                       = output_folder;
  if (!utils::folder_exists(output_folder))
    utils::mkdir(output_folder.c_str());

  if (paired_output)
    this->init_all_file_handles();
  else
    this->init_all_file_handles_single_end();
}


// return and open file handles
void
FileHandler::init_all_file_handles()
{
  string barcode, sample_name;

  for (auto it = Barcode_file_map->begin(); it != Barcode_file_map->end(); it++) {
    barcode     = it->first;
    sample_name = it->second;
    string                                        mate1_path = Output_folder + "/" + sample_name +
                                                               "_R1.fastq.gz";
    string                                        mate2_path = Output_folder + "/" + sample_name +
                                                               "_R2.fastq.gz";
    ZipFastqWriter                                *w1       = new ZipFastqWriter(mate1_path);
    ZipFastqWriter                                *w2       = new ZipFastqWriter(mate2_path);
    std::pair<ZipFastqWriter *, ZipFastqWriter *> *filepair = new std::pair<ZipFastqWriter *,
                                                                            ZipFastqWriter *>(w1,
                                                                                              w2);
    Fastq_handler[barcode]                = filepair;
    Map_filehandle_sample_name[filepair]  = sample_name;
  }
}


void
FileHandler::init_all_file_handles_single_end()
{
  string          barcode, sample_name;
  ZipFastqWriter  *dummy_writer = NULL;

  for (auto it = Barcode_file_map->begin(); it != Barcode_file_map->end(); it++) {
    barcode     = it->first;
    sample_name = it->second;
    string                                        mate1_path = Output_folder + "/" + sample_name +
                                                               ".fastq.gz";
    ZipFastqWriter                                *w1       = new ZipFastqWriter(mate1_path);
    std::pair<ZipFastqWriter *, ZipFastqWriter *> *filepair = new std::pair<ZipFastqWriter *,
                                                                            ZipFastqWriter *>(w1,
                                                                                              dummy_writer);
    Fastq_handler[barcode]                = filepair;
    Map_filehandle_sample_name[filepair]  = sample_name;
  }
}


std::pair<ZipFastqWriter *, ZipFastqWriter *> *
FileHandler::get_file_handles(string &barcode)
{
  auto it = Fastq_handler.find(barcode);

  if (it != Fastq_handler.end())
    return it->second;
  else
    return Fastq_handler["undetermined"];

  return NULL;
}


string
FileHandler::get_sample_name(std::pair<ZipFastqWriter *, ZipFastqWriter *> *file_pair_key)
{
  auto it = Map_filehandle_sample_name.find(file_pair_key);

  if (it != Map_filehandle_sample_name.end())
    return it->second;

  return "";
}


FileHandler::~FileHandler()
{
  // clear file handles
  for (auto it = Fastq_handler.begin(); it != Fastq_handler.end(); it++) {
    if (it->second->first)
      delete it->second->first;

    if (it->second->second)
      delete it->second->second;

    delete it->second;
  }
}
