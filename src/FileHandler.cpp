/*
 * FileHandler.cpp
 *
 *  Created on: Aug 15, 2020
 *      Author: quin
 */

#include "FileHandler.h"
#include <unordered_map>
#include <string>
#include <cmath>
#include <sys/stat.h>

#include <sstream>
#include <sys/stat.h>

// for windows mkdir
#ifdef _WIN32
#include <direct.h>
#endif

namespace utils
{
    /**
     * Checks if a folder exists
     * @param foldername path to the folder to check.
     * @return true if the folder exists, false otherwise.
     */
    bool folder_exists(std::string foldername)
    {
        struct stat st;
        stat(foldername.c_str(), &st);
        return st.st_mode & S_IFDIR;
    }

    /**
     * Portable wrapper for mkdir. Internally used by mkdir()
     * @param[in] path the full path of the directory to create.
     * @return zero on success, otherwise -1.
     */
    int _mkdir(const char *path)
    {
    #ifdef _WIN32
        return ::_mkdir(path);
    #else
    #if _POSIX_C_SOURCE
        return ::mkdir(path, 0755);
    #else
        return ::mkdir(path, 0755); // not sure if this works on mac
    #endif
    #endif
    }

    /**
     * Recursive, portable wrapper for mkdir.
     * @param[in] path the full path of the directory to create.
     * @return zero on success, otherwise -1.
     */
    int mkdir(const char *path)
    {
        std::string current_level = "";
        std::string level;
        std::stringstream ss(path);

        // split path using slash as a separator
        while (std::getline(ss, level, '/'))
        {
            current_level += level; // append folder to the current level

            // create current level
            if (!folder_exists(current_level) && _mkdir(current_level.c_str()) != 0)
                return -1;

            current_level += "/"; // don't forget to append a slash
        }

        return 0;
    }
}

FileHandler::FileHandler(unordered_map<string,string>& barcode_file_map, string output_folder, size_t memory) {
  // TODO Auto-generated constructor stub
	Barcode_file_map = &barcode_file_map;
	Output_folder = output_folder;
	//Fastq_handler = new unordered_map<ZipFastqWriter,ZipFastqWriter>();
	int buffer_per_file = int(memory / (barcode_file_map.size()) * 2);

	if (!utils::folder_exists(output_folder))
		utils::mkdir(output_folder.c_str());
	init_all_file_handles();
}

// return and open file handles
void FileHandler::init_all_file_handles(){
	string barcode, sample_name;
	for(auto it = Barcode_file_map->begin(); it != Barcode_file_map->end(); it++){
		barcode = it->first;
		sample_name = it->second;
		string mate1_path = Output_folder +"/"+ sample_name + "_R1.fastq.gz";
		string mate2_path = Output_folder +"/"+ sample_name + "_R2.fastq.gz";
		std::pair<ZipFastqWriter*,ZipFastqWriter*> filepair = {new ZipFastqWriter(mate1_path),new ZipFastqWriter(mate2_path)};
		Fastq_handler[barcode] = filepair;
	}
}

std::pair<ZipFastqWriter*,ZipFastqWriter*> FileHandler::get_file_handles(string &barcode){
	return Fastq_handler[barcode];
}

/*
 class FileHandler(ExitStack):

    def __init__(self, barcode_file_map, output_folder, memory=2 ** 30):
        self.barcode_file_map = barcode_file_map
        self.output_folder = pathlib.Path(output_folder)
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        self.buffer_per_file = int(memory // (len(barcode_file_map) * 2))
        super().__init__()
        self.fastq_handler = defaultdict(list)
        self.barcode_file_map["undetermined"] = "undetermined"

    def __enter__(self):
        for barcode, sample_name in self.barcode_file_map.items():
            mate1_path = self.output_folder / "{}{}".format(sample_name, "_R1.fastq.gz")
            mate2_path = self.output_folder / "{}{}".format(sample_name, "_R2.fastq.gz")
            self.fastq_handler[barcode] = self._open_gz_file_handles(
                mate1_path, mate2_path
            )
        return self.fastq_handler

    def _open_gz_file_handles(self, mate1_path, mate2_path):
        out = []
        for mate_path in (mate1_path, mate2_path):
            gz_out = io.BufferedWriter(gzip.open(mate_path, mode='wb', compresslevel=4),
                                       buffer_size=self.buffer_per_file)
            out.append(gz_out)
            log.debug("File handler opened for mate: %s", mate_path)
        return out
 */


FileHandler::~FileHandler() {
  // TODO Auto-generated destructor stub
	// clear file handles
	for(auto it = Fastq_handler.begin(); it != Fastq_handler.end(); it++){
		delete it->second.first;
		delete it->second.second;
	}
}

