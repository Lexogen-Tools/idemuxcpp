#include <thread>
#include <math.h>
#include "config.h"
#include "PairedReader.h"
#include "FastqReader.h"
#include "BoostZipReader.h"
#if HAVE_LIBBAMTOOLS
#include "BamReader.h"
#endif

PairedReader::PairedReader(string fastqgz_1,
                           string fastqgz_2,
                           double queue_buffer_gb) : Queue_buffer_bytes(max(size_t(1),
                                                                            size_t (queue_buffer_gb
                                                                                    * pow(1024,
                                                                                          3)))),
  IsInterleaved(false)
{
  if (fastqgz_1 == "")
    throw runtime_error("Error: read 1 file name is empty!");

  if (fastqgz_2 == "")
    throw runtime_error("Error: read 1 file name is empty!");

  string  suffix_r1     = "";
  string  suffix_r2     = "";
  size_t  ext_index_r1  = fastqgz_1.rfind(".");

  if (ext_index_r1 != fastqgz_1.npos)
    suffix_r1 = fastqgz_1.substr(ext_index_r1);

  size_t  ext_index_r2 = fastqgz_2.rfind(".");

  if (ext_index_r2 != fastqgz_2.npos)
    suffix_r2 = fastqgz_2.substr(ext_index_r2);

  if (suffix_r1.compare(".gz") == 0)
    // use zip reader
    Reader1 = new BoostZipReader(fastqgz_1);

  string error_bam =
    string("Error: not interleaved bam files are not supported with a read 2 file!") + \
    string(
      " If it is 1 bam file that contains reads in PE mode, please use only read 1 file and the paired flag!");

  if (suffix_r1.compare(".bam") == 0)
    throw runtime_error(error_bam);

  if (suffix_r1.compare(".gz") != 0 && suffix_r1.compare(".bam") != 0)
    // use fastq reader
    Reader1 = new FastqReader(fastqgz_1);

  if (suffix_r2.compare(".gz") == 0)
    // use zip reader
    Reader2 = new BoostZipReader(fastqgz_2);

  if (suffix_r2.compare(".bam") == 0)
    throw runtime_error(error_bam);

  if (suffix_r2.compare(".gz") != 0 && suffix_r2.compare(".bam") != 0)
    // use fastq reader
    Reader2 = new FastqReader(fastqgz_2);
}


/**
 * Constructor for single end mode or interleaved paired mode.
 * @param single_read_file - read 1 file or interleaved paired end reads file.
 * @param queue_buffer_gb - buffer for reads queue.
 * @param is_interleaved - if false, assume single end mode, else paired end mode (the only PE option for .bam files is interleaved).
 */
PairedReader::PairedReader(string single_read_file,
                           double queue_buffer_gb,
                           bool   is_interleaved) : Queue_buffer_bytes(max(size_t(1),
                                                                           size_t (queue_buffer_gb *
                                                                                   pow(1024, 3))))
{
  string  suffix    = "";
  size_t  ext_index = single_read_file.rfind(".");

  if (ext_index != single_read_file.npos)
    suffix = single_read_file.substr(ext_index);

  if (suffix.compare(".gz") == 0)
    // use zip reader
    Reader1 = new BoostZipReader(single_read_file);

  if (suffix.compare(".bam") == 0) {
#if HAVE_LIBBAMTOOLS
    Reader1 = new BamReader(single_read_file);
#else
    throw runtime_error("Error: bam files are not supported with this compilation!");
#endif
  }

  if (suffix.compare(".gz") != 0 && suffix.compare(".bam") != 0)
    // use fastq reader
    Reader1 = new FastqReader(single_read_file);

  IsInterleaved = is_interleaved;
  Reader2       = NULL;
}


void
PairedReader::fill_reads_interleaved(IFastqReader                                   *reader,
                                     size_t                                         max_size,
                                     std::vector<std::pair<fq_read *, fq_read *> >  *pairs)
{
  size_t bytes_read = 0;

  while ((Queue_buffer_bytes > 0 && (bytes_read < Queue_buffer_bytes)) ||
         (max_size > 0 && (pairs->size() < max_size))) {
    std::pair<fq_read *, fq_read *> pair = reader->next_read_pair();
    if (pair.first != NULL && pair.second != NULL) {
      pairs->push_back(pair);
      if (Queue_buffer_bytes > 0) {
        bytes_read  += pair.first->get_size();
        bytes_read  += pair.second->get_size();
      }
    } else {
      break;
    }
  }
}


std::vector<std::pair<fq_read *, fq_read *> > *
PairedReader::next_reads(size_t max_size)
{
  std::vector<std::pair<fq_read *, fq_read *> > *pairs = new std::vector<std::pair<fq_read *,
                                                                                   fq_read *> >();

  //fq_read* pairs = new fq_read[max_size*2+1];
  if (IsInterleaved and Reader2 == NULL) {
    fill_reads_interleaved(Reader1, max_size, pairs);
  } else {
    fq_read *read1, *read2;
    size_t  bytes_read = 0;
    while ((Queue_buffer_bytes > 0 && (bytes_read < Queue_buffer_bytes)) ||
           (max_size > 0 && (pairs->size() < max_size))) {
      read1 = Reader1->next_read();
      if (Reader2 != NULL)
        read2 = Reader2->next_read();
      else
        read2 = NULL; // dummy read in single end mode.

      if (read1 != NULL) {
        pairs->push_back(std::pair<fq_read *, fq_read *>(read1, read2));
        if (Queue_buffer_bytes > 0) {
          bytes_read += read1->get_size();
          if (read2 != NULL)
            bytes_read += read2->get_size();
        }
      } else {
        delete read1;
        delete read2;
        break;
      }
    }
  }

  return pairs;
}


void
fill_reads(IFastqReader *reader,
           size_t       max_size,
           fq_read      ***reads)
{
  for (size_t i = 0; i < max_size; i++) {
    fq_read *tmp = reader->next_read();
    if (tmp != NULL) {
      (*reads)[i] = tmp;
    } else {
      (*reads)[i] = NULL;
      delete tmp;
      break;
    }
  }
}


std::vector<std::pair<fq_read *, fq_read *> > *
PairedReader::next_reads2(size_t  max_size,
                          int     reading_threads)
{
  std::vector<std::pair<fq_read *, fq_read *> > *pairs = new std::vector<std::pair<fq_read *,
                                                                                   fq_read *> >();

  if (IsInterleaved and Reader2 == NULL) {
    fill_reads_interleaved(Reader1, max_size, pairs);
  } else {
    if (reading_threads == 2) {
      fq_read     **reads1  = new fq_read *[max_size]{ 0 };
      fq_read     **reads2  = new fq_read *[max_size]{ 0 };

      std::thread to1(fill_reads, Reader1, max_size, &reads1);
      if (Reader2 != NULL) {
        std::thread to2(fill_reads, Reader2, max_size, &reads2);
        to2.join();
      }

      to1.join();

      for (size_t i = 0; i < max_size; i++) {
        if (reads1[i] != NULL)
          pairs->push_back(std::pair<fq_read *, fq_read *>(reads1[i], reads2[i]));
        else
          break;
      }
      delete[] reads1;
      delete[] reads2;
    } else {
      // use 1 thread
      pairs = this->next_reads(max_size);
    }
  }

  return pairs;
}


PairedReader::~PairedReader()
{
  delete Reader1;
  delete Reader2;
}
