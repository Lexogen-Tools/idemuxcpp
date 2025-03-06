#include <string>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include "Writer.h"
#include "Parser.h"

using namespace std;

Writer::Writer()
{
}


/**
 * Writes a counter dictionary as tsv file with the name error_correction_stats.tsv.
 *
 * @param counter(dict): A dictionary barcodes <sample name, #corrected reads>.
 * @param output_dir (string): The path the file should be written to.
 *
 */
void
Writer::write_summary(std::unordered_map<string, size_t> &counter,
                      string                              output_dir)
{
  string    output_file = string("demultipexing_stats.tsv");

  if (output_dir.length() > 0)
    output_file = output_dir + PATH_SEP + string("demultipexing_stats.tsv");

  ofstream  csvfile(output_file, std::ofstream::out);
  string    csv_header = "sample_name\twritten_reads";

  csvfile << csv_header << std::endl;
  string    sample_filename;
  size_t    counts;

  vector<string> sorted_file_names;
  for (auto it = counter.begin(); it != counter.end(); it++) {
      sorted_file_names.push_back(it->first);
  }
  sort(sorted_file_names.begin(), sorted_file_names.end());

  for (auto it = sorted_file_names.begin(); it != sorted_file_names.end(); it++) {
    sample_filename = *it;
    counts          = counter[sample_filename];
    csvfile << sample_filename << "\t" << counts << std::endl;
  }
  fprintf(stdout, "Run complete! Summary statistics saved to %s\n",
          output_file.c_str());
}


struct {
  bool
  operator()(const pair<string, size_t> & a,
             const pair<string, size_t> & b) const
  {
    if (a.second == b.second)
      return a.first.compare(b.first);

    return a.second < b.second;
  }
} compare_barcode_count;


Writer::~Writer()
{
}
