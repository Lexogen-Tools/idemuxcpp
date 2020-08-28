#include <string>
#include <stdio.h>
#include <unordered_map>
#include "Writer.h"

using namespace std;

Writer::Writer() {
  // TODO Auto-generated constructor stub

}

/**
 * Writes a counter dictionary as tsv file with the name error_correction_stats.tsv.
 *
 * @param counter(dict): A dictionary barcodes <sample name, #corrected reads>.
 * @param output_dir (string): The path the file should be written to.
 *
 */
void Writer::write_summary(std::unordered_map<string, size_t> counter, string output_dir){
      string output_file = output_dir + "/demultipexing_stats.tsv";
      char open_mode = 'w';
      FILE *csvfile = fopen(output_file.c_str(), &open_mode);
      string csv_header = "sample\twritten_reads";
      fwrite(csv_header.c_str(),sizeof(char),csv_header.length(),csvfile);
      string sample_filename;
      size_t counts;
      for(auto it = counter.begin(); it != counter.end(); it++){
        sample_filename = it->first;
        counts = it->second;
        fprintf(csvfile, string(sample_filename + "\t%ud").c_str(),counts);
      }
      fprintf(stdout, "Run complete! Summary statistics saved to %s", output_file.c_str());
      fclose(csvfile);
}

Writer::~Writer() {
  // TODO Auto-generated destructor stub
}

