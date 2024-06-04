#ifndef PARSER_H_
#define PARSER_H_

#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <set>
#include "Barcode.h"
#include "helper.h"
#include "FastqReader.h"
#include "BoostZipReader.h"
#include "PairedReader.h"

#ifdef _WIN32
#include "direct.h"
#define PATH_SEP '\\'
#define GETCWD _getcwd
#define CHDIR _chdir
#else
#include "unistd.h"
#define PATH_SEP '/'
#define GETCWD getcwd
#define CHDIR chdir
#endif
#include <sys/stat.h>

using namespace std;


struct i1_info {
  // read 1 or read 2
  int read_index;
  // start index (1-based)
  int start_index;
  // i1 barcode is within [start_index, end_index].
  int end_index;
};

class Parser {
public:
Parser();
string
remove_bom(string s);


std::vector<std::vector<std::string> >
readCSV(std::istream &in);


unordered_map<string, string> *
parse_sample_sheet(string                           sample_sheet,
                   bool                             i5_rc,
                   bool                             i7_rc,
                   bool                             auto_detect,
                   vector<Barcode *> &              barcodes_out,
                   unordered_map<string, i1_info> & i7_i5_i1_info_map,
                   string                           relative_exepath,
                   string                           correction_maps_path,
                   bool                             demux_only,
                   int                              default_i1_read,
                   int                              default_i1_start,
                   bool                             single_end_mode = false);


string
reverse_complement(string sequence);


bool
has_valid_barcode_combinations(Barcode &i7,
                               Barcode &i5,
                               Barcode &i1);


void
fastq_lines_to_reads(string fastq_lines);


void
peek_into_fastq_files(PairedReader &                  get_pe_fastq,
                      bool                            has_i7,
                      bool                            has_i5,
                      bool                            has_i1,
                      vector<int> &                   i7_length,
                      vector<int> &                   i5_length,
                      unordered_map<string, i1_info> &i7_i5_i1_info_map,
                      size_t                          max_length_i7 = SIZE_MAX,
                      size_t                          max_length_i5 = SIZE_MAX);


void
peek_into_fastq_file(IFastqReader                     *reader,
                     bool                             has_i7,
                     bool                             has_i5,
                     bool                             has_i1,
                     vector<int> &                    i7_length,
                     vector<int> &                    i5_length,
                     unordered_map<string, i1_info> & i7_i5_i1_info_map,
                     size_t                           max_length_i7 = SIZE_MAX,
                     size_t                           max_length_i5 = SIZE_MAX);


void
check_mate_pair(std::pair<fq_read *, fq_read *> mate_pair,
                bool                            has_i7,
                bool                            has_i5,
                bool                            has_i1,
                vector<int> &                   i7_length,
                vector<int> &                   i5_length,
                unordered_map<string, i1_info> &i7_i5_i1_info_map,
                size_t                          max_length_i7 = SIZE_MAX,
                size_t                          max_length_i5 = SIZE_MAX);


void
check_mate2_length(fq_read  *mate2,
                   int      i1_start,
                   int      i1_end);


void
check_fastq_headers(std::pair<fq_read *, fq_read *> mate_pair,
                    bool                            has_i7,
                    bool                            has_i5,
                    vector<int> &                   i7_length,
                    vector<int> &                   i5_length,
                    size_t                          max_length_i7 = SIZE_MAX,
                    size_t                          max_length_i5 = SIZE_MAX);


void
check_fastq_header(fq_read      *mate,
                   bool         has_i7,
                   bool         has_i5,
                   vector<int> &i7_length,
                   vector<int> &i5_length,
                   size_t       max_length_i7 = SIZE_MAX,
                   size_t       max_length_i5 = SIZE_MAX);


static std::pair<string, string>
parse_indices(const string &input,
              size_t        max_length_i7 = SIZE_MAX,
              size_t        max_length_i5 = SIZE_MAX)
{
  size_t  index_colon = input.find_last_of(':');
  size_t  index_plus1 = input.find_last_of('+');
  string  code_i7     = "";
  string  code_i5     = "";

  if (index_colon != std::string::npos)
    code_i7 = input.substr(index_colon + 1, min(index_plus1 - index_colon - 1, max_length_i7));

  if (index_plus1 != std::string::npos)
    code_i5 = input.substr(index_plus1 + 1, min(input.length() - index_plus1 - 1, max_length_i5));

  std::pair<string, string> bcs_mate1 = std::pair<string, string>(
    code_i7,
    code_i5);

  return bcs_mate1;
}


static string
replace_indices(const string &input,
                const string &i7,
                const string &i5)
{
  size_t  index_colon = input.find_last_of(':');
  string  result      = "";

  if (index_colon != std::string::npos) {
    result.append(input.substr(0, index_colon + 1));
    result.append(i7);
    if (i7.length() > 0 and i5.length() > 0)
      result.append("+");

    result.append(i5);
  } else {
    result = input;
  }

  return result;
}


static
bool
PathExists(const std::string &s)
{
  struct stat buffer;

  return stat(s.c_str(), &buffer) == 0;
}


/*
 * load error correction map from misc folder.
 */
static
void
get_map_from_resource(string                        filepath,
                      unordered_map<string, string> *mapping,
                      set<string>                   *barcodes)
{
  printf("Loading error correction map from %s\n", filepath.c_str());

  ifstream  dataFile;

  dataFile.open(filepath);
  fprintf(stdout, "checking file for valid barcodes: %s\n", filepath.c_str());
  int       tabidx;
  int       return_idx; //for windows files
  string    v1, v2;

  if (dataFile.fail()) {
    string message = string_format(
      "Error: failed to open the correction map file! Please copy it "
      "to the following location: %s\n", filepath.c_str());
    throw(runtime_error(message));
  }

  while (!dataFile.eof()) {
    std::string temp = "";
    while (std::getline(dataFile, temp)) {
      tabidx      = temp.find('\t');
      return_idx  = temp.find('\r');
      if (!return_idx)
        return_idx = temp.length();

      if (tabidx >= 0 && tabidx < (int)temp.length()) {
        v1  = temp.substr(0, tabidx);
        v2  = temp.substr(tabidx + 1, return_idx - tabidx - 1);
        //std::cout << "row" << v1 << ", " << v2 << std::endl;
        mapping->insert({ v1, v2 });
        barcodes->insert(v2);
      }

      temp = "";
    }
  }
}


virtual
~Parser();
private:
string
list_to_string(vector<int> &list);
};

#endif /* PARSER_H_ */
