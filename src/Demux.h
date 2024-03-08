#ifndef SRC_DEMUX_H_
#define SRC_DEMUX_H_

#include "config.h"
#include <zlib.h>
#include <iostream>
#include <queue>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <limits>
#include <cmath>
#include <chrono>

#include "Parser.h"
#include "FileHandler.h"
#include "PairedReader.h"
#include "FastqReader.h"
#if HAVE_LIBBAMTOOLS
#include "BamReader.h"
#endif
#include "BoostZipReader.h"
#include "ZipFastqWriter.h"
#include "Writer.h"
#include "Correction_Counter.h"

#include <future>
#include <thread>

using namespace std;

typedef std::chrono::high_resolution_clock Clock;


string
extract_i1(fq_read  *r,
           size_t   i1_start,
           size_t   i1_end)
{
  return r->Sequence.substr(i1_start, i1_end - i1_start);
}


/**
 * remove the i1 interval from the sequence and the quality code and append it to the index.
 */
void
cut_out_i1(fq_read  *r,
           size_t   i1_start,
           size_t   i1_end)
{
  string i1_bc = r->Sequence.substr(i1_start, i1_end - i1_start);

  r->Seq_ID   = string(r->Seq_ID).append("+").append(i1_bc);
  r->Sequence = string(r->Sequence).substr(0, i1_start).append(r->Sequence.substr(i1_end));
  //r->Plus_ID = string(r->Plus_ID);
  r->QualityCode = string(r->QualityCode).substr(0, i1_start).append(r->QualityCode.substr(i1_end));
}


/**
 * process mate pair
 * ignore read 2 if it is NULL. Then it is single end mode.
 * Returns
 * -------
 * barcodes: ()
 */
string
process_mate_pair(std::pair<fq_read *, fq_read *> & mate_pair,
                  unordered_set<string>             *i7_wanted,
                  unordered_set<string>             *i5_wanted,
                  unordered_set<string>             *i1_wanted,
                  unordered_map<string, string>     *map_i7,
                  unordered_map<string, string>     *map_i5,
                  unordered_map<string, string>     *map_i1,
                  unordered_map<string, i1_info> &  i7_i5_i1_info_map,
                  std::pair<fq_read *, fq_read *> & mate_pair_out,
                  Correction_Counter                *counted_corrections_per_index,
                  size_t                            max_length_i7 = SIZE_MAX,
                  size_t                            max_length_i5 = SIZE_MAX)
{
  string                    fastq_header(mate_pair.first->Seq_ID);
  std::pair<string, string> barcodes = Parser::parse_indices(fastq_header,
                                                             max_length_i7,
                                                             max_length_i5);
  // if not barcodes:
  string                    i7_bc = "";
  string                    i5_bc = "";

  // when there are 2 barcodes in the fastq header the orientation is i7,i5
  string                    tmp_i7_bc = barcodes.first;
  string                    tmp_i5_bc = barcodes.second;

  if (i7_wanted->size() == 0)
    i7_bc = "";

  if (i5_wanted->size() == 0)
    i5_bc = "";

  bool  corrected_i7, corrected_i5, corrected_i1;

  corrected_i7 = corrected_i5 = corrected_i1 = false;

  auto  i7_bc_it = map_i7->find(tmp_i7_bc);

  if (i7_bc_it != map_i7->end()) {
    i7_bc = i7_bc_it->second;
    if (counted_corrections_per_index != NULL && i7_bc.compare(tmp_i7_bc) != 0)
      corrected_i7 = true;
  }

  auto i5_bc_it = map_i5->find(tmp_i5_bc);

  if (i5_bc_it != map_i5->end()) {
    i5_bc = i5_bc_it->second;
    if (counted_corrections_per_index != NULL && i5_bc.compare(tmp_i5_bc) != 0)
      corrected_i5 = true;
  }

  string  i1_bc = "";

  fq_read *r1c  = new fq_read(*mate_pair.first);
  fq_read *r2c  = NULL;

  if (mate_pair.second != NULL) // if PE mode.
    r2c = new fq_read(*mate_pair.second);

  /* correct i7 and i5 indices (write the original indices,
   * that could have a different length after parsing because of
   * the --restrict-barcode-length option)
   */
  string new_seq_id = Parser::replace_indices(r1c->Seq_ID, tmp_i7_bc, tmp_i5_bc);

  r1c->Seq_ID = new_seq_id;
  if (r2c != NULL) {
    new_seq_id  = Parser::replace_indices(r2c->Seq_ID, tmp_i7_bc, tmp_i5_bc);
    r2c->Seq_ID = new_seq_id;
  }

  if (i1_wanted->size() > 0) {
    string  i7_i5_bc    = i7_bc + "\n" + i5_bc;
    auto    it_i1_info  = i7_i5_i1_info_map.find(i7_i5_bc);
    if (it_i1_info != i7_i5_i1_info_map.end()) {
      size_t  i1_start  = it_i1_info->second.start_index;
      size_t  i1_end    = it_i1_info->second.end_index;
      // correct only i1 with length > 0
      if (i1_start < i1_end) {
        if (it_i1_info->second.read_index == 1)
          i1_bc = extract_i1(r1c, i1_start, i1_end);
        else
          if (r2c != NULL)
            i1_bc = extract_i1(r2c, i1_start, i1_end);

        auto iti1c = map_i1->find(i1_bc);
        if (iti1c != map_i1->end()) {
          string  _i1_corrected = iti1c->second;
          //if _i1_corrected in i1_wanted:
          auto    it1 = i1_wanted->find(_i1_corrected);
          if (it1 != i1_wanted->end()) {
            if (it_i1_info->second.read_index == 1) {
              cut_out_i1(r1c, i1_start, i1_end);
              if (r2c != NULL)
                r2c->Seq_ID.append("+").append(i1_bc);
            } else {
              if (r2c != NULL) {
                cut_out_i1(r2c, i1_start, i1_end);
                r1c->Seq_ID.append("+").append(i1_bc);
              }
            }

            if (counted_corrections_per_index != NULL && i1_bc.compare(_i1_corrected) != 0)
              corrected_i1 = true;

            i1_bc = _i1_corrected;
          }
        }
      }
    }
  }

  mate_pair_out.first   = r1c;
  mate_pair_out.second  = r2c;
  string barcodes_key = string("").append(i7_bc).append("\n").append(i5_bc).append("\n").append(
    i1_bc);

  if (counted_corrections_per_index != NULL)
    counted_corrections_per_index->count_correction(i7_bc,
                                                    i5_bc,
                                                    i1_bc,
                                                    corrected_i7,
                                                    corrected_i5,
                                                    corrected_i1);

  return barcodes_key;
}


size_t
get_max_length(unordered_set<string> &barcodes)
{
  size_t  max       = 0;
  size_t  bc_length = 0;

  for (auto it = barcodes.begin(); it != barcodes.end(); it++) {
    bc_length = (*it).length();
    if (bc_length > max)
      max = bc_length;
  }
  return max;
}


typedef struct demux_thread_parameter_ {
  std::vector<std::pair<fq_read *, fq_read *> >             *pe_reads; //input
  unordered_set<string>                                     *i7_wanted;
  unordered_set<string>                                     *i5_wanted;
  unordered_set<string>                                     *i1_wanted;
  unordered_map<string, string>                             *map_i7;
  unordered_map<string, string>                             *map_i5;
  unordered_map<string, string>                             *map_i1;
  unordered_map<string, i1_info>                            *i7_i5_i1_info_map;
  size_t                                                    max_length_i7;
  size_t                                                    max_length_i5;
  FileHandler                                               *file_handler;
  unordered_map<std::pair<ZipFastqWriter *, ZipFastqWriter *> *,
                vector<std::pair<fq_read *, fq_read *> > *> *map_pairs;                     // output
  Correction_Counter                                        *counted_corrections_per_index; // output
} demux_thread_parameter;


int
demux_list_reads(demux_thread_parameter *tp)
{
  tp->map_pairs = new unordered_map<std::pair<ZipFastqWriter *, ZipFastqWriter *> *,
                                    vector<std::pair<fq_read *, fq_read *> > *>();
  size_t                                        max_reads = tp->pe_reads->size();
  std::pair<fq_read *, fq_read *>               mate_pair;
  std::pair<ZipFastqWriter *, ZipFastqWriter *> *pair_writer;
  string                                        s_barcodes;

  for (size_t i = 0; i < max_reads; i++) {
    mate_pair = tp->pe_reads->at(i);
    std::pair<fq_read *, fq_read *> processed_mates;
    s_barcodes = process_mate_pair(mate_pair,
                                   tp->i7_wanted,
                                   tp->i5_wanted,
                                   tp->i1_wanted,
                                   tp->map_i7,
                                   tp->map_i5,
                                   tp->map_i1,
                                   *tp->i7_i5_i1_info_map,
                                   processed_mates,
                                   tp->counted_corrections_per_index,
                                   tp->max_length_i7,
                                   tp->max_length_i5);
    pair_writer = tp->file_handler->get_file_handles(s_barcodes);

    //output:
    if (!(*tp->map_pairs)[pair_writer])
      (*tp->map_pairs)[pair_writer] = new vector<std::pair<fq_read *, fq_read *> >();

    (*tp->map_pairs)[pair_writer]->push_back(processed_mates);
  }
  return 0;
}


int
compress_and_delete(string          *string_to_write,
                    vector<uint8_t> *compressed_memory)
{
  ZipFastqWriter::compress_memory((uint8_t *)string_to_write->c_str(),
                                  string_to_write->length(),
                                  compressed_memory,
                                  4);
  delete string_to_write;
  return 0;
}


typedef std::pair<ZipFastqWriter *, size_t> compression_input;

void
demux(unordered_map<string, string>   *barcode_sample_map,
      vector<Barcode *> &             barcodes,
      string                          read1_file,
      string                          read2_file,
      unordered_map<string, i1_info> &i7_i5_i1_info_map,
      string                          output_dir,
      Parser &                        parser,
      size_t                          queue_size,
      size_t                          reading_threads,
      size_t                          compressing_threads,
      size_t                          processing_threads,
      string                          barcode_corrections_file,
      bool                            skip_check,
      bool                            restrict_barcode_length,
      double                          queue_buffer_gb,
      size_t                          gzip_block_size_bytes,
      bool                            single_end_mode,
      bool                            verbose)
{
  // load the maps that will be used for error correction. As the tool does not allow
  // different length we only need to load the used length
  Barcode *i7, *i5, *i1;

  i7  = barcodes[0];
  i5  = barcodes[1];
  i1  = barcodes[2];

  unordered_set<string>         *i7_wanted  = i7->used_codes();
  unordered_set<string>         *i5_wanted  = i5->used_codes();
  unordered_set<string>         *i1_wanted  = i1->used_codes();
  unordered_map<string, string> *map_i7     = &i7->Correction_map;
  unordered_map<string, string> *map_i5     = &i5->Correction_map;
  unordered_map<string, string> *map_i1     = &i1->Correction_map;

  // if None is in *_wanted no barcode has been specified
  bool                          has_i7  = !i7->empty();
  bool                          has_i5  = !i5->empty();
  bool                          has_i1  = !i1->empty();

  size_t                        max_length_i7 = SIZE_MAX;
  size_t                        max_length_i5 = SIZE_MAX;

  if (restrict_barcode_length) {
    max_length_i7 = get_max_length(*i7_wanted);
    max_length_i5 = get_max_length(*i5_wanted);
  }

  bool  is_interleaved  = (read2_file == "" && !single_end_mode);
  bool  paired_output   = (read2_file != "") || is_interleaved;

  // before doing any processing check if the fastq file is okay.
  if (!skip_check) {
    PairedReader *pr;
    if (read2_file == "")
      //if read1_file is not interleaved, it could be paiered end mode.
      pr = new PairedReader(read1_file, queue_buffer_gb, is_interleaved);
    else
      pr = new PairedReader(read1_file, read2_file, queue_buffer_gb);

    parser.peek_into_fastq_files(*pr,
                                 has_i7,
                                 has_i5,
                                 has_i1,
                                 i7->Lengths,
                                 i5->Lengths,
                                 i7_i5_i1_info_map,
                                 max_length_i7,
                                 max_length_i5);
    delete pr;
  }

  Correction_Counter *counted_corrections_per_index = NULL;

  if (barcode_corrections_file.compare("") != 0)
    counted_corrections_per_index = new Correction_Counter(i7_wanted, i5_wanted, i1_wanted);

  //read_counter = Counter()
  std::cout << "Starting demultiplexing" << std::endl;
  PairedReader *pr;

  if (read2_file == "")
    //if read1_file is not interleaved, it could be paiered end mode.
    pr = new PairedReader(read1_file, queue_buffer_gb, is_interleaved);
  else
    pr = new PairedReader(read1_file, read2_file, queue_buffer_gb);

  // first we need to open the output files the reads should get sorted into
  FileHandler                                                 *file_handler
    = new FileHandler(*barcode_sample_map, output_dir, paired_output);

  // then we iterate over all the paired end reads
  std::vector<std::pair<fq_read *,
                        fq_read *> >                          *pe_reads = NULL;
  std::unordered_map<std::string, size_t>                     read_counter;

  unordered_map<demux_thread_parameter *, std::future<int> >  demux_threads;

  unordered_map<std::pair<compression_input, vector<uint8_t> *> *,
                std::future<int> >                            compress_threads;

  unordered_map<ZipFastqWriter *,
                string *>                                     uncompressed_read_string_buffer;
  unordered_map<ZipFastqWriter *,
                vector<string *> >                            uncompressed_read_string_chunks;

  size_t
                                                              compressed_block_id = 1;
  size_t
                                                              last_compressed_block_id =
    compressed_block_id - 1;
  /* id for compressed block.
   * Needed to maintain Read1 and Read2 reads order in final files.
   * Sequential write happens after compression pool.
   * It is an edge case, but we must take care of it.
   */

  double  reads_demuxed         = 0;
  double  reads_read            = 0;
  bool    finished_all_reading  = false;
  bool    finish_while          = false;
  auto    t1                    = std::chrono::high_resolution_clock::now();

  while (!finish_while) {
    size_t uncompressed_chunks = 0;
    for (auto it = uncompressed_read_string_chunks.begin();
         it != uncompressed_read_string_chunks.end(); it++)
      uncompressed_chunks += it->second.size();

    if (pe_reads == NULL && compress_threads.size() < compressing_threads &&
        uncompressed_chunks < compressing_threads) {
      // only read when writing threads are not overloaded (because processing is fast and would accumulate to much memory for results).
      pe_reads    = pr->next_reads2(queue_size, reading_threads);
      reads_read  += pe_reads->size();
      if (pe_reads->size() == 0) {
        delete pe_reads;
        pe_reads              = NULL;
        finished_all_reading  = true;
      }
    }

    vector<demux_thread_parameter *> dt_to_erase;
    for (auto it_dt = demux_threads.begin(); it_dt != demux_threads.end(); it_dt++) {
      bool is_done = true;
      try {
        if (it_dt->second.valid())
          is_done = it_dt->second.wait_for(std::chrono::nanoseconds(0)) ==
                    std::future_status::ready;
      } catch (std::future_error &e) {
        is_done = false;
        std::cerr << "Error: " << e.what() << std::endl;
      }
      if (is_done) {
        // collect thread results in queue for demuxed but uncompressed reads.
        demux_thread_parameter *tp = it_dt->first;
        dt_to_erase.push_back(tp);
        if (tp != NULL) {
          if (tp->map_pairs != NULL) {
            for (auto it_mapped = tp->map_pairs->begin(); it_mapped != tp->map_pairs->end();
                 it_mapped++) {
              //results.push_back(std::pair(it_mapped->first, it_mapped->second));
              vector<std::pair<fq_read *,
                               fq_read *> > *reads = it_mapped->second;
              std::pair<ZipFastqWriter *,
                        ZipFastqWriter *>   *writer_pair          = it_mapped->first;
              string                        *uncompressed_reads1  = NULL;
              string                        *uncompressed_reads2  = NULL;
              if (uncompressed_read_string_buffer.find(writer_pair->first) ==
                  uncompressed_read_string_buffer.end()) {
                uncompressed_reads1                                 = new string("");
                uncompressed_read_string_buffer[writer_pair->first] = uncompressed_reads1;
              }

              if (writer_pair->second != NULL) {
                if (uncompressed_read_string_buffer.find(writer_pair->second) ==
                    uncompressed_read_string_buffer.end()) {
                  uncompressed_reads2                                   = new string("");
                  uncompressed_read_string_buffer[writer_pair->second]  = uncompressed_reads2;
                }
              }

              uncompressed_reads1 = uncompressed_read_string_buffer[writer_pair->first];
              if (writer_pair->second != NULL)
                uncompressed_reads2 = uncompressed_read_string_buffer[writer_pair->second];

              for (auto it = reads->begin(); it != reads->end(); it++) {
                if (uncompressed_reads1->size() < gzip_block_size_bytes) {
                  uncompressed_reads1->append(it->first->to_string());
                } else {
                  uncompressed_read_string_chunks[writer_pair->first].push_back(uncompressed_reads1);
                  uncompressed_reads1 = new string(
                    it->first->to_string());
                  uncompressed_read_string_buffer[writer_pair->first] = uncompressed_reads1;
                }

                if (it->second != NULL) {
                  if (uncompressed_reads2->size() < gzip_block_size_bytes) {
                    uncompressed_reads2->append(it->second->to_string());
                  } else {
                    uncompressed_read_string_chunks[writer_pair->second].push_back(
                      uncompressed_reads2);
                    uncompressed_reads2 = new string(
                      it->second->to_string());
                    uncompressed_read_string_buffer[writer_pair->second] = uncompressed_reads2;
                  }
                }
              }

              // count and erase
              size_t  number_reads  = reads->size();
              string  sample_name   = file_handler->get_sample_name(writer_pair);
              auto    it_rc         = read_counter.find(sample_name);
              if (it_rc == read_counter.end())
                read_counter[sample_name] = number_reads;
              else
                it_rc->second += number_reads;

              reads_demuxed += number_reads;

              // erase processed reads.
              for (auto it = reads->begin(); it != reads->end(); it++) {
                delete it->first;
                delete it->second;
              }
              delete reads;
            }
            delete tp->map_pairs;
            // merge corrections
            if (barcode_corrections_file.compare("") != 0) {
              *counted_corrections_per_index += *tp->counted_corrections_per_index;
              delete tp->counted_corrections_per_index;
            }
          }

          // clean up input
          for (auto it = tp->pe_reads->begin(); it != tp->pe_reads->end(); it++) {
            delete it->first;
            delete it->second;
          }
          delete tp->pe_reads;
        }
      }
    }
    // clear threads that are done.
    for (size_t i = 0; i < dt_to_erase.size(); i++) {
      demux_threads.erase(dt_to_erase[i]);
      delete dt_to_erase[i];
    }
    dt_to_erase.clear();

    // start new thread if reads are there and a thread in the pool is empty.
    if (pe_reads != NULL && demux_threads.size() < processing_threads) {
      demux_thread_parameter *tp = new demux_thread_parameter();
      tp->pe_reads                      = pe_reads;
      tp->i7_wanted                     = i7_wanted;
      tp->i5_wanted                     = i5_wanted;
      tp->i1_wanted                     = i1_wanted;
      tp->map_i7                        = map_i7;
      tp->map_i5                        = map_i5;
      tp->map_i1                        = map_i1;
      tp->max_length_i7                 = max_length_i7;
      tp->max_length_i5                 = max_length_i5;
      tp->i7_i5_i1_info_map             = &i7_i5_i1_info_map;
      tp->file_handler                  = file_handler; // separate filehandler and reset pointers for writing.
      tp->counted_corrections_per_index = NULL;         // separate map & merge afterwards.
      if (barcode_corrections_file.compare("") != 0)
        tp->counted_corrections_per_index = new Correction_Counter(i7_wanted, i5_wanted, i1_wanted);

      demux_threads[tp] = std::async(std::launch::async, demux_list_reads, tp);
      pe_reads          = NULL;
    }

    // take demuxed results from queue and start separate threads
    vector<std::pair<compression_input, vector<uint8_t> *> *> ct_to_erase;
    for (auto it_ct = compress_threads.begin(); it_ct != compress_threads.end(); it_ct++) {
      bool is_done = true;
      try {
        if (it_ct->second.valid()) {
          is_done = it_ct->second.wait_for(
            std::chrono::nanoseconds(0))
                    == std::future_status::ready;
        }
      } catch (std::future_error &e) {
        // not ready
        is_done = false;
      }
      if (is_done) {
        // compressing done.
        // check if it is the next block, otherwise leave it there and continue.
        bool    is_next_block = false;
        size_t  block_id      = it_ct->first->first.second;
        if (block_id == last_compressed_block_id + 1)
          is_next_block = true;

        if (is_next_block) {
          // -> write results.
          std::pair<compression_input, vector<uint8_t> *> *cr     = it_ct->first;
          compression_input                               ci      = cr->first;
          ZipFastqWriter                                  *writer = ci.first;
          last_compressed_block_id = ci.second;
          vector<uint8_t>                                 *compressed_data = cr->second;
          writer->write(compressed_data);
          delete compressed_data;
          ct_to_erase.push_back(cr);
          auto                                            t2 =
            std::chrono::high_resolution_clock::now();
          double                                          seconds =
            std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
          if (verbose) {
            string output_pattern = string(
              "\rreads: %10.4g, reads demuxed: %10.4g, seconds: %10.4g, reads per second: %10.4g, ")
                                    +
                                    string(
              "uncompressed chunks: %lu, demux threads: %lu, compress threads: %lu") +
                                    std::string(to_string(UINT64_MAX).length() * 2 -
                                                to_string(uncompressed_chunks).length() -
                                                to_string(seconds).length(), ' ').c_str();

            printf(output_pattern.c_str(),
                   reads_read,
                   reads_demuxed,
                   seconds,
                   reads_demuxed / seconds,
                   uncompressed_chunks,
                   demux_threads.size(),
                   compress_threads.size());
          } else {
            printf(
              "\rreads: %10.4g, reads demuxed: %10.4g, seconds: %10.4g, reads per second: %10.4g",
              reads_read,
              reads_demuxed,
              seconds,
              reads_demuxed / seconds);
          }

          fflush(stdout);
        }
      }
    }
    // erase finished threads
    for (size_t i = 0; i < ct_to_erase.size(); i++) {
      compress_threads.erase(ct_to_erase[i]);
      delete ct_to_erase[i];
    }
    ct_to_erase.clear();

    //check if we compress last chunk of demuxed data
    if (finished_all_reading && demux_threads.size() == 0 &&
        uncompressed_read_string_chunks.size() == 0 && uncompressed_read_string_buffer.size() > 0) {
      for (auto it_uc = uncompressed_read_string_buffer.begin();
           it_uc != uncompressed_read_string_buffer.end(); it_uc++)
        uncompressed_read_string_chunks[it_uc->first].push_back(it_uc->second);
      uncompressed_read_string_buffer.clear();
    }

    // for each writer pair start a thread until thread pool is full.
    while (uncompressed_read_string_chunks.size() > 0 &&
           compress_threads.size() < compressing_threads) {
      auto it_first_res = uncompressed_read_string_chunks.begin();
      // write compressed data.
      if (it_first_res != uncompressed_read_string_chunks.end()) {
        ZipFastqWriter    *writer = it_first_res->first;
        vector<string *>  *reads  = &it_first_res->second;
        if (reads->size() > 0) {
          vector<uint8_t>               *compressed_memory =
            new vector<uint8_t>();
          compression_input             compression_in =
          { writer, compressed_block_id };
          std::pair<compression_input,
                    vector<uint8_t> *>  *compressed_result = new std::pair<compression_input,
                                                                           vector<uint8_t> *>(
            { compression_in, compressed_memory });
          string                        *uc_data = reads->at(0);
          compress_threads[compressed_result] = std::async(std::launch::async,
                                                           compress_and_delete,
                                                           uc_data,
                                                           compressed_memory);
          compressed_block_id++;
          reads->erase(reads->begin());
        } else {
          uncompressed_read_string_chunks.erase(writer);
        }
      }
    }
    finish_while = demux_threads.size() == 0 && compress_threads.size() == 0 &&
                   finished_all_reading && uncompressed_read_string_buffer.size() == 0 &&
                   uncompressed_read_string_chunks.size() == 0;
  }
  std::cout << std::endl;
  if (barcode_corrections_file.compare("") != 0) {
    counted_corrections_per_index->write_correction(barcode_corrections_file, *barcode_sample_map);
    delete counted_corrections_per_index;
  }

  Writer summary_writer;

  summary_writer.write_summary(read_counter, output_dir);

  delete pr;
  delete file_handler;
  delete i7_wanted;
  delete i5_wanted;
  delete i1_wanted;
}


#endif /* SRC_DEMUX_H_ */
