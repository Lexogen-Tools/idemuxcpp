/*
 * ZipFastqWriter.h
 *
 *  Created on: Aug 26, 2020
 *      Author: gentzian
 */

#ifndef ZIPFASTQWRITER_H_
#define ZIPFASTQWRITER_H_

#include <string>
#include "ZipFastqReader.h"

using namespace std;

class ZipFastqWriter {
public:
	ZipFastqWriter(string outfile);
	ZipFastqWriter(const ZipFastqWriter&) = delete; // keine Kopie
	ZipFastqWriter& operator=(const ZipFastqWriter&) = delete; // keine Zuweisung
	int write(const char* data, size_t len);
	void write_read(fq_read* r);
	void write_read_list(fq_read** reads);
	virtual ~ZipFastqWriter();
private:
	string OutputFile;
	gzFile GZ_file_handle;
};

#endif /* ZIPFASTQWRITER_H_ */
