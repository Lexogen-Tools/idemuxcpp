#ifndef WRITER_H_
#define WRITER_H_

#include <unordered_map>

class Writer {
public:
  Writer();
  void write_summary(std::unordered_map<std::string, size_t> &counter, std::string output_dir);
  virtual ~Writer();
};

#endif /* WRITER_H_ */

