#ifndef WRITER_H_
#define WRITER_H_

#include <unordered_map>
#include <string>

using namespace std;

struct correction_counter {
  unordered_map<string, size_t> count_i1;
  unordered_map<string, size_t> count_i5;
  unordered_map<string, size_t> count_i7;
};

class Writer {
public:
Writer();
void
write_summary(unordered_map<string, size_t> & counter,
              string                          output_dir);

virtual
~Writer();
};

#endif /* WRITER_H_ */
