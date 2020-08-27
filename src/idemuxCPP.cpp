//============================================================================
// Name        : idemuxCPP.cpp
// Author      : GE
// Version     :
// Copyright   : (C) GE 2020
// Description : Demultiplex fastq.gz reads according to their indices.
//============================================================================

// compile with: g++ -std=gnu++11 -o idemuxCPP idemuxCPP.cpp -lz

#include <zlib.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <cstring>
#include <algorithm>
#include <unordered_map>
#include <limits>
#include <cmath>

#include<idemuxCPP_cmdl.h>
#include "Parser.h"
#include "FileHandler.h"
#include "PairedReader.h"
#include "ZipFastqReader.h"
#include "ZipFastqWriter.h"

using namespace std;


/**
 * """process mate pair

    Returns
    -------
    barcodes: ()
    """
 */
string process_mate_pair(std::pair<fq_read*,fq_read*>& mate_pair,
		unordered_set<string>* i7_wanted, unordered_set<string>* i5_wanted, unordered_set<string>* i1_wanted,
                      bool has_i7,
					  unordered_map<string,string> &map_i7, unordered_map<string,string> &map_i5, unordered_map<string,string> &map_i1,
                      int i1_start, int i1_end, std::pair<fq_read*,fq_read*>& mate_pair_out){


    string fastq_header(mate_pair.first->Seq_ID);
    std::pair<string,string> barcodes = Parser::parse_indices(fastq_header);
    // if not barcodes:
    string i7_bc, i5_bc;

    // when there are 2 barcodes in the fastq header the orientation is i7,i5
    //if len(barcodes) == 2:
    //    i7_bc, i5_bc = barcodes
    i7_bc = barcodes.first;
    i5_bc = barcodes.second;

    //elif len(barcodes) == 1:
    //    i7_bc, i5_bc = (barcodes[0], None) if has_i7 else (None, barcodes[0])
    if(!has_i7){
    	i5_bc = i7_bc;
    	i7_bc = "";
    }

    i7_bc = map_i7.at(i7_bc);
    i5_bc = map_i5.at(i5_bc);

    string i1_bc = "";
    auto it7 = i7_wanted->find(i7_bc);
    auto it5 = i5_wanted->find(i5_bc);

    //if i7_bc in i7_wanted and i5_bc in i5_wanted:
    if (it7 != i7_wanted->end() and it5 != i5_wanted->end()){
    	string i1_bc = string(mate_pair.second->Sequence).substr(i1_start,i1_end);
        //i1_bc = mate_pair[1][1][i1_start:i1_end]
    	string _i1_corrected = map_i1.at(i1_bc);
    	//_i1_corrected = map_i1.get(i1_bc)

        //if _i1_corrected in i1_wanted:
    	auto it1 = i1_wanted->find(_i1_corrected);
    	if(it1 != i1_wanted->end()){
            i1_bc = _i1_corrected;

            //(m1_hdr, m1_seq, m1_opt, m1_qcs), (m2_hdr, m2_seq, m2_opt, m2_qcs) = mate_pair
            fq_read* r1c = new fq_read();
            fq_read* r2c = new fq_read();
            r1c->Seq_ID = string(mate_pair.first->Seq_ID)+ "+" + i1_bc;
    	    //m1_hdr = f"{m1_hdr[:-1]}+{i1_bc}\n"
            r1c->Sequence = string(mate_pair.first->Sequence);
            r1c->Plus_ID= string(mate_pair.first->Plus_ID);
            r1c->QualityCode= string(mate_pair.first->QualityCode);

            r2c->Seq_ID = string(mate_pair.first->Seq_ID)+ "+" + i1_bc;
			r2c->Sequence = string(mate_pair.first->Sequence).substr(0,i1_start) + string(mate_pair.first->Sequence).substr(i1_end,mate_pair.first->Sequence.length());
			r2c->Plus_ID= string(mate_pair.first->Plus_ID);
			r2c->QualityCode= string(mate_pair.first->QualityCode).substr(0,i1_start) + string(mate_pair.first->QualityCode).substr(i1_end,mate_pair.first->QualityCode.length());
    	    //m2_hdr = f"{m2_hdr[:-1]}+{i1_bc}\n"
    	    //m2_seq = f"{m2_seq[:i1_start]}{m2_seq[i1_end:]}"
    	    //m2_qcs = f"{m2_qcs[:i1_start]}{m2_qcs[i1_end:]}"

			mate_pair_out.first = r1c;
			mate_pair_out.second = r2c;
    	}
    }
    string barcodes_key = i7_bc +'\n'+ i5_bc +'\n'+ i1_bc;
    return barcodes_key;
}

/*

def process_mate_pair(mate_pair,
                      i7_wanted, i5_wanted, i1_wanted,
                      has_i7,
                      map_i7, map_i5, map_i1,
                      i1_start, i1_end):
    """process mate pair

    Returns
    -------
    barcodes: ()
    """

    fastq_header = mate_pair[0][0]
    _, _, _barcodes = fastq_header.rpartition(":")
    barcodes = _barcodes[:-1].split("+")
    # if not barcodes:
    i7_bc = i5_bc = None

    # when there are 2 barcodes in the fastq header the orientation is i7,i5
    if len(barcodes) == 2:
        i7_bc, i5_bc = barcodes

    elif len(barcodes) == 1:
        i7_bc, i5_bc = (barcodes[0], None) if has_i7 else (None, barcodes[0])

    i7_bc = map_i7.get(i7_bc)
    i5_bc = map_i5.get(i5_bc)

    i1_bc = None

    if i7_bc in i7_wanted and i5_bc in i5_wanted:

        i1_bc = mate_pair[1][1][i1_start:i1_end]
        _i1_corrected = map_i1.get(i1_bc)

        if _i1_corrected in i1_wanted:
            i1_bc = _i1_corrected
            (m1_hdr, m1_seq, m1_opt, m1_qcs), (m2_hdr, m2_seq, m2_opt, m2_qcs) = mate_pair

            mate_1 = (
                f"{m1_hdr[:-1]}+{_i1_corrected}\n"
                f"{m1_seq[:i1_start]}{m1_seq[i1_end:]}"
                f"{m1_opt}"
                f"{m1_qcs[:i1_start]}{m1_qcs[i1_end:]}"
            )

            mate_2 = (
                f"{m2_hdr[:-1]}+{_i1_corrected}\n"
                f"{m2_seq[:i1_start]}{m2_seq[i1_end:]}"
                f"{m2_opt}"
                f"{m2_qcs[:i1_start]}{m2_qcs[i1_end:]}"
            )

            mate_pair = (mate_1, mate_2)
    return (i7_bc, i5_bc, i1_bc), mate_pair


def demux_paired_end(barcode_sample_map, barcodes, read1, read2, i1_start, output_dir,
                     **kwargs):
    # TODO: add documentation
    # TODO: add logging
    # load the maps that will be used for error correction. As the tool does not allow
    # different length we only need to load the used length
    i7, i5, i1 = barcodes
    i7_wanted = i7.used_codes
    i5_wanted = i5.used_codes
    i1_wanted = i1.used_codes
    map_i7 = i7.correction_map
    map_i5 = i5.correction_map
    map_i1 = i1.correction_map
    i1_end = i1_start + i1.length

    # if None is in *_wanted no barcode has been specified
    has_i7 = not i7.empty
    has_i5 = not i5.empty
    has_i1 = not i1.empty

    # before doing any processing check if the fastq file is okay.
    peek_into_fastq_files(read1, read2, has_i7, has_i5, has_i1, i7.length,
                          i5.length, i1_start, i1_end)
    read_counter = Counter()
    log.info("Staring demultiplexing")
    # first we need to open the output files the reads should get sorted into
    with FileHandler(barcode_sample_map, output_dir) as file_handler:
        # then we iterate over all the paired end reads
        with get_pe_fastq(read1, read2) as pe_reads:
            for mate_pair in tqdm(pe_reads):
                # here we do the error correction and get obtain the i1 barcode if present
                barcodes, processed_mates = process_mate_pair(mate_pair,
                                                              i7_wanted,
                                                              i5_wanted,
                                                              i1_wanted,
                                                              has_i7,
                                                              map_i7,
                                                              map_i5,
                                                              map_i1,
                                                              i1_start,
                                                              i1_end)
                # When a barcode combination is unknown/not specified in the sample sheet
                # it will get  sorted into file containing only reads with unknown
                # barcodes.
                fq_out1, fq_out2 = file_handler.get(barcodes,
                                                    file_handler["undetermined"])
                fq_out1.write(processed_mates[0].encode())
                fq_out2.write(processed_mates[1].encode())
    write_summary(read_counter, output_dir)
 */
void demux_paired_end(unordered_map<string,string>* barcode_sample_map, vector<Barcode*>& barcodes, string read1, string read2, int i1_start, string output_dir, Parser& parser){
    // TODO: add documentation
    // TODO: add logging
    // load the maps that will be used for error correction. As the tool does not allow
    // different length we only need to load the used length
    Barcode *i7, *i5, *i1;
    i7 = barcodes[0];
    i5 = barcodes[1];
    i1 = barcodes[2];

    unordered_set<string>* i7_wanted = i7->used_codes();
    unordered_set<string>* i5_wanted = i5->used_codes();
    unordered_set<string>* i1_wanted = i1->used_codes();
    unordered_map<string,string> map_i7 = i7->correction_map;
    unordered_map<string,string> map_i5 = i5->correction_map;
    unordered_map<string,string> map_i1 = i1->correction_map;
    int i1_end = i1_start + i1->length;

    // if None is in *_wanted no barcode has been specified
    bool has_i7 = !i7->empty();
    bool has_i5 = !i5->empty();
    bool has_i1 = !i1->empty();

    // before doing any processing check if the fastq file is okay.
    parser.peek_into_fastq_files(read1, read2, has_i7, has_i5, has_i1, i7->length,
                          i5->length, i1_start, i1_end);

    //read_counter = Counter()
    printf("Staring demultiplexing");
    // first we need to open the output files the reads should get sorted into
    FileHandler *file_handler = new FileHandler(*barcode_sample_map, output_dir,(size_t)pow(2,30));
    //with FileHandler(barcode_sample_map, output_dir) as file_handler:
        // then we iterate over all the paired end reads
    PairedReader pr(read1,read2);
    std::vector<std::pair<fq_read*,fq_read*>>* pe_reads;
    bool reads_available = true;
    while(reads_available){
    	pe_reads = pr.next_reads(1);
    	if(pe_reads == NULL){
    		reads_available = false;
    		break;
    	}
    	for(auto it = pe_reads->begin(); it != pe_reads->end(); it++){
    		std::pair<fq_read*,fq_read*> mate_pair = *it;
    		std::pair<fq_read*,fq_read*> processed_mates;
    		string barcodes = process_mate_pair(mate_pair,
									  i7_wanted,
									  i5_wanted,
									  i1_wanted,
									  has_i7,
									  map_i7,
									  map_i5,
									  map_i1,
									  i1_start,
									  i1_end, processed_mates);
    		std::pair<ZipFastqWriter*,ZipFastqWriter*> pair_writer = file_handler->get_file_handles(barcodes);
    		//processed_mates
    		pair_writer.first->write_read(processed_mates.first);
    		pair_writer.first->write_read(processed_mates.second);
    	}
    }
/*
        with get_pe_fastq(read1, read2) as pe_reads:
            for mate_pair in tqdm(pe_reads):
                // here we do the error correction and get obtain the i1 barcode if present
                barcodes, processed_mates = process_mate_pair(mate_pair,
                                                              i7_wanted,
                                                              i5_wanted,
                                                              i1_wanted,
                                                              has_i7,
                                                              map_i7,
                                                              map_i5,
                                                              map_i1,
                                                              i1_start,
                                                              i1_end)
                // When a barcode combination is unknown/not specified in the sample sheet
                // it will get  sorted into file containing only reads with unknown
                // barcodes.
                fq_out1, fq_out2 = file_handler->  get(barcodes,
                                                    file_handler["undetermined"])
                fq_out1.write(processed_mates[0].encode())
                fq_out2.write(processed_mates[1].encode())
                */
    //write_summary(read_counter, output_dir);

}

void usage() {
  printf(
      "usage:\nidemuxCPP -t <transcript size> -f <fragment size> -r <read size> --mode <1 for spp, 2 for fcp, 3 for rcp>\n");
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {
  string read1_file = "";
  string read2_file = "";
  string outputdirectory = "";
  string sample_sheet_file = "";
  bool i5_rc = false;
  int i1_start = 10;
  bool verbose = false;

  struct idemuxCPP_args_info args_info;
  // check if there is any unparseable argument
  if (idemuxCPP_cmdline_parser(argc, argv, &args_info) != 0)
    exit(EXIT_FAILURE);

  if (args_info.r1_given) {
    read1_file = args_info.r1_arg;
  } else {
    fprintf(stderr, "Error: please enter a read1.fastq.gz file!\n");
    exit(EXIT_FAILURE);
  }
  if (args_info.r2_given) {
    read1_file = args_info.r2_arg;
  } else {
    fprintf(stderr, "Error: please enter a read2.fastq.gz file!\n");
    exit(EXIT_FAILURE);
  }
  if (args_info.out_given) {
    outputdirectory = args_info.out_arg;
  }
  if (args_info.sample_sheet_given) {
    sample_sheet_file = args_info.sample_sheet_arg;
  }
  i5_rc = args_info.i5_rc_flag;
  i1_start = args_info.i1_start_arg;
  verbose = args_info.verbose_flag;


  Parser p;
  vector<Barcode*> barcodes;
  unordered_map<string,string>* barcode_sample_map = p.parse_sample_sheet(sample_sheet_file, i5_rc, barcodes);

  // do things.
  demux_paired_end(barcode_sample_map, barcodes, read1_file, read2_file, i1_start, outputdirectory,p);

  idemuxCPP_cmdline_parser_free(&args_info);
  return EXIT_SUCCESS;
}

