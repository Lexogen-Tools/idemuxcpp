#ifndef ZIPFASTQWRITER_H_
#define ZIPFASTQWRITER_H_

#include <string>
#include <vector>
#include <zlib.h>
#include "FastqReader.h"

using namespace std;

class ZipFastqWriter {
public:
ZipFastqWriter(string outfile);
ZipFastqWriter(const ZipFastqWriter&) = delete;   //no copy
ZipFastqWriter&
operator=(const ZipFastqWriter&) = delete;        // no assignment


int
write(std::vector<uint8_t> *out_data);


virtual
~ZipFastqWriter();
inline string
get_output_name()
{
  return OutputFile;
}


static int
compress_memory(uint8_t               *in_data,
                size_t                in_data_size,
                std::vector<uint8_t>  *out_data,
                int                   compression_level);


private:
string OutputFile;
FILE *File_handle;
};

#endif /* ZIPFASTQWRITER_H_ */
