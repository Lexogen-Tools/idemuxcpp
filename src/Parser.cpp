/*
 * Parser.cpp
 *
 *  Created on: Aug 15, 2020
 *      Author: quin
 */

#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include "Parser.h"
#include "Barcode.h"
#include "PairedReader.h"

#include <istream>
#include <fstream>
#include <string>
#include <vector>

enum class CSVState {
    UnquotedField,
    QuotedField,
    QuotedQuote
};

std::vector<std::string> readCSVRow(const std::string &row) {
    CSVState state = CSVState::UnquotedField;
    std::vector<std::string> fields {""};
    size_t i = 0; // index of the current field
    for (char c : row) {
        switch (state) {
            case CSVState::UnquotedField:
                switch (c) {
                    case ',': // end of field
                              fields.push_back(""); i++;
                              break;
                    case '"': state = CSVState::QuotedField;
                              break;
                    default:  fields[i].push_back(c);
                              break; }
                break;
            case CSVState::QuotedField:
                switch (c) {
                    case '"': state = CSVState::QuotedQuote;
                              break;
                    default:  fields[i].push_back(c);
                              break; }
                break;
            case CSVState::QuotedQuote:
                switch (c) {
                    case ',': // , after closing quote
                              fields.push_back(""); i++;
                              state = CSVState::UnquotedField;
                              break;
                    case '"': // "" -> "
                              fields[i].push_back('"');
                              state = CSVState::QuotedField;
                              break;
                    default:  // end of quote
                              state = CSVState::UnquotedField;
                              break; }
                break;
        }
    }
    return fields;
}

/// Read CSV file, Excel dialect. Accept "quoted fields ""with quotes"""
std::vector<std::vector<std::string>> readCSV(std::istream &in) {
    std::vector<std::vector<std::string>> table;
    std::string row;
    while (!in.eof()) {
        std::getline(in, row);
        if (in.bad() || in.fail()) {
            break;
        }
        auto fields = readCSVRow(row);
        table.push_back(fields);
    }
    return table;
}

Parser::Parser() {
  // TODO Auto-generated constructor stub

}

/**
 * Function to parse the sample_sheet.csv.

      This function takes the path to a idemux sample sheet reads the data and does some
      sanity checks to prevent downstream problems. The csv needs to consist out of 4
      columns specifing a sample name and the respective barcodes. If a sample is only
      i7, i5 or i1 barcoded missing barcodes should be indicated by am empty field.
      In general the sample sheet should be formatted like this and have a header included:

      sample_name,i7,i5,i1
      <sample_name>,<i7_barcode or empty>,<i7_barcode or empty>,<i1_barcode or empty>

      Example:
      These are allowed
          triple indexing (all barcodes specified)
          sample_name,i7,i5,i1
          sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
          sample_1,AAAATCCCAGTT,CCCCTAAACGTT,AAAATCCCAGTT

          dual indexing(either i7+i5, or i7+i1 or i5+i1)
          sample_name,i7,i5,i1
          sample_0,AAAACATGCGTT,,AAAACATGCGTT
          sample_1,AAAATCCCAGTT,,AAAATCCCAGTT

          single indexing(either i7, i5 or i5)
          sample_name,i7,i5,i1
          sample_0,,,AAAACATGCGTT
          sample_1,,,AAAATCCCAGTT

          mixed i1 indexing(some samples have dual or triple indexing. Works only with i1!)
          sample_name,i7,i5,i1
          sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
          sample_1,AAAATCCCAGTT,CCCCTAAACGTT,

          different length for either all i7 (10 nt) or i5 (12 nt)
          sample_name,i7,i5,i1
          sample_0,AAAACATGCG,CCCCACTGAGTT,AAAACATGCGTT
          sample_1,AAAATCCCAG,CCCCTAAACGTT,AAAATCCCAGTT


      These are not allowed and will call sys.exit(1):
          mixed i7/i5 indexing
          sample_name,i7,i5,i1
          sample_0,AAAACATGCGTT,,AAAACATGCGTT
          sample_1,,CCCCTAAACGTT,AAAATCCCAGTT

          duplicate/non-unique barcodes
          sample_name,i7,i5,i1
          sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
          sample_1,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT

          duplicate sample names
          sample_name,i7,i5,i1
          sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
          sample_0,AAAATCCCAGTT,CCCCTAAACGTT,AAAATCCCAGTT

          different length for one barcode
          sample_name,i7,i5,i1
          sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATGCGTT
          sample_1,AAAATCCCAG,CCCCTAAACGTT,AAAATCCCAGTT

          an i1 barcode shorter than 12 nt
          sample_name,i7,i5,i1
          sample_0,AAAACATGCGTT,CCCCACTGAGTT,AAAACATG
          sample_1,AAAATCCCAG,CCCCTAAACGTT,AAAATCCC

      Arguments:
          sample_sheet (str): Path to a sample_sheet.cvs. See examples for formatting

      Returns:
          barcode_sample_map (dict): A dict mapping a tuple of barcodes to sample names
              (i7, i5, i1): sample_name
          used_barcodes (tuple): A tuple of 3 sets, each containing the i7, i5 and i1
              barcodes specified in the sample sheet
          barcode_lengths (dict): A dict mapping the index name to a list barcode
              lengths.  {"i7": list(int), "i5": list(int), "i1": list(int)}

      Except:
          ValueError: Will initiate sys.exit(1)
 */
unordered_map<string,string>* Parser::parse_sample_sheet(string sample_sheet, bool i5_rc, vector<Barcode*>& barcodes_out) {
      // we use these to keep track which lengths are beeing used for each barcode
      unordered_set<int> i7_lengths, i5_lengths, i1_lengths;

      // we use these to keep track of which single match to which sample(s)
      unordered_map<string,std::vector<string>> i7_barcodes;
      unordered_map<string,std::vector<string>> i5_barcodes;
      unordered_map<string,std::vector<string>> i1_barcodes;

      /* this is for keeping track which unique barcode combination is associated with a sample */
      unordered_map<string,string>* barcode_sample_map = new unordered_map<string,string>();
      // we use this for checking for duplicated sample names
      unordered_map<string,size_t> sample_count;
      // expected file header
      //TODO: use a set, such that the columns dont have to be in this order.
      vector<string> sample_sheet_header = {"sample_name", "i7", "i5", "i1"};
      //char open_mode = 'r';

      ifstream in_file_sample_sheet(sample_sheet.c_str());
      // we register a csv dialect to do whitespace trimming for us
      //csv.register_dialect('strip', skipinitialspace=True)
      size_t  lines_read = 0;
      //reader = csv.DictReader(sample_sheet, restval=None, dialect='strip')
      std::vector<std::vector<std::string>> csv_lines = readCSV(in_file_sample_sheet);

      if(csv_lines.size() == 0){
    	  fprintf(stderr, "Error: the sample sheet is empty!");
		exit(EXIT_FAILURE);
      }
      //check header:
      string s1,s2;
      //std::cout << "sizeof " << sample_sheet_header.size() << std::endl;
      for(int i = 0; i< sample_sheet_header.size(); i++){
        s1 = csv_lines[0][i];
        s2 = sample_sheet_header[i];
        if (s1.compare(s2) != 0){
          fprintf(stderr, "Incorrect sample sheet header. Expected header: %s\n"
                                                     "Observed header: %s",s1.c_str(), s2.c_str());
        }
      }
      // check file
      for(size_t i = 1; i < csv_lines.size(); i++){
        std::vector<std::string> row = csv_lines[i];

        lines_read += 1;
        string sample_name = row[0];

        // we use this to keep track of duplicate sample names
        sample_count[sample_name] += 1;
        string tmp_i7 = row[1];
        string i7_bc = tmp_i7; // has i7.
        if (tmp_i7.compare("") == 0){
          i7_bc = "None";
        }

        // i5 can be sequenced as reverse complement, translate if needed
        string tmp_i5 = row[2];
        if (i5_rc){
          tmp_i5 = reverse_complement(row[2]);
        }
        string i5_bc = tmp_i5;
        if (tmp_i5.compare("") == 0){
          i5_bc = "None";
        }

        string tmp_i1 = row[3];
        string i1_bc = tmp_i1; // has i1.
        if (tmp_i1.compare("") == 0){
          i1_bc = "None";
        }
        //std::cout << sample_name << tmp_i7 <<  tmp_i5 << tmp_i1 << std::endl;
        // add barcodes and sample_names to dict so we can do some value checking
        // later as not all combinations should be allowed
        i7_barcodes[i7_bc].push_back(sample_name);
        i5_barcodes[i5_bc].push_back(sample_name);
        i1_barcodes[i1_bc].push_back(sample_name);
        string barcodes = i7_bc +'\n'+ i5_bc +'\n'+ i1_bc;
        // barcode combinations have to be unique. otherwise we don't know to which
        // sample they belong
        auto it_bc = barcode_sample_map->find(barcodes);
        if (it_bc != barcode_sample_map->end()){
          string same_barcode = barcode_sample_map->at(barcodes);
          fprintf(stderr, "Duplicate barcode combination detected. Sample: %s "
                                           "Barcodes: %s\n Already observed for: %s\n Each "
                                           "barcode combination has to be unique.", sample_name,
                                                                                      barcodes,
                                                                                     same_barcode);
          exit(EXIT_FAILURE);
        }

        (*barcode_sample_map)[barcodes] = sample_name;
        // barcodes of different length are a problem as 8, 10, 12 nucleotide
        // barcodes are a subset of each other. Therefore we only allow one
        // length per barcode type
        if(i7_bc.compare("None") != 0)
          i7_lengths.insert(i7_bc.length());
        if(i5_bc.compare("None") != 0)
                  i7_lengths.insert(i5_bc.length());
        if(i5_bc.compare("None") != 0)
                  i7_lengths.insert(i5_bc.length());
      }
         in_file_sample_sheet.close();


         Barcode* i7 = new Barcode("i7", i7_barcodes);
          Barcode* i5 = new Barcode("i5", i5_barcodes, i5_rc);
          Barcode* i1 = new Barcode("i1", i1_barcodes);

          /*
          barcodes_out.clear();
          barcodes_out.push_back(i7);
          barcodes_out.push_back(i5);
          barcodes_out.push_back(i1);
	*/

          // test if the supplied barcode combinations are valid

         bool is_valid = has_valid_barcode_combinations(*i7, *i5, *i1);
          // sample names have to be unique as they determine the outfile name. Otherwise
          // we get problems when we try to write reads belonging to different barcode
          // combinations to one file.
          //TODO: correct this
          bool duplicated_sample_names = false; //[k for k, v in sample_count.items() if v > 1]
          string duplicated_sample_names_str = "";
          for(auto it = sample_count.begin(); it != sample_count.end(); it++){
        	  if (it->second > 1){
        		  duplicated_sample_names = true;
        		  duplicated_sample_names_str += it->first + ", ";
        	  }
          }
          if (duplicated_sample_names){

              string error_msg = "The sample sheet contains duplicate sample names. Sample "\
                           "names have to be unique!\n Sample names: "+\
						   duplicated_sample_names_str;
              fprintf(stderr,error_msg.c_str());
              exit(EXIT_FAILURE);
          }

      // we have made it until here until raising an exception. That means the sample sheet
      // information should be okay and we can return required data from the sample sheet.
      //TODO:
      //barcodes = [load_correction_map(bc) for bc in barcodes]
          Barcode* bc7 = load_correction_map(*i7);
          Barcode* bc5 = load_correction_map(*i5);
          Barcode* bc1 = load_correction_map(*i1);
          barcodes_out.push_back(bc7);
          barcodes_out.push_back(bc5);
          barcodes_out.push_back(bc1);
      //return barcode_sample_map, barcodes
          return barcode_sample_map;
}

/**
 * Function that returns the reverse complement of DNA sequence. Accepts A,C,T,G,N
    as input bases.

    Args:
        sequence (str): A DNA string that should be translated to its reverse complement

    Returns:
        str: Reverse complement of the input. Returns None when
            sequence is None

    Raises:
        ValueError: Is raised when the input string contains other letters than A,C,T,G,N
*/
string Parser::reverse_complement(string sequence) {
    // if sequence is None, no need to do any work
    if (sequence == "")
        return "";

    // complements of each base
    unordered_map<char,char> base_complements;
    base_complements['A'] = 'T';
    base_complements['C'] = 'G';
    base_complements['G'] = 'C';
    base_complements['T'] = 'A';
    base_complements['N'] = 'N';
    try
    {
        // get the reverse complement
        string rc_sequence = "";
        for(int i = 0; i < sequence.length(); i++){
        	auto it_comp = base_complements.find(sequence[i]);
        	if (it_comp != base_complements.end())
        		rc_sequence += it_comp->second;
        }
        //= [base_complements[nt] for nt in sequence[::-1]];
        return rc_sequence; //"".join(rc_sequence);
    } catch (const std::exception& e) {
        // we only get a key error if sequence contains letter that are not covered
        // above
    	string invalid_bases = "";
		for(int i = 0; i < sequence.length(); i++){
			auto it_comp = base_complements.find(sequence[i]);
			if (it_comp == base_complements.end())
				invalid_bases += sequence[i];
		}
        string error_message = "The following barcode sequence from the sample sheet "\
                         "contains bases that can't be mapped to their reverse "\
                         "complement.\nBarcodes: %s\nBases %s";
        fprintf(stderr, error_message.c_str(), sequence.c_str(), invalid_bases.c_str());
        exit(EXIT_FAILURE);
    };
    return "";
}

/*
 * Function that checks if the provided barcodes allow unique sample identification
    for the different usecases.

    Allowed barcoding variants are:
    - all contain an i7 and i5
    - all contain no i7 all an i5
    - all contain an i7 no an i5
    - all contain no i7, no i5, all an i1

    Forbidden barcoding variants are:
    - some have an i7
    - some have an i5
    - no i7 no i5, no i1

    Args:
        i7 (Barcode): Barcode dataclass of i7 barcodes
        i5 (Barcode): Barcode dataclass of i5 barcodes
        i1 (Barcode): Barcode dataclass of i1 barcodes

    Returns:
        is_valid (bool): True when valid barcode combinations are supplied (see above)

    Raises:
        ValueError (err): When supplied an invalid barcode combination (see above)

 */
bool Parser::has_valid_barcode_combinations(Barcode& i7, Barcode& i5, Barcode& i1) {
   /*
    allowed cases, we dont actually need these return values, but these make the logic
    much more obvious and human readable.
    demultiplexing on i7 and i5 (and maybe i1). I1 is optional when i7 and i5 are
    already specified so we dont really need to check for sparsity
    */
    if (i7.full() && i5.full())
        return true;
    // demultiplexing on i7 and/or i1
    if (i7.full() && i5.empty())
        return true;
    // demultiplexing on i5 and/or i1
    if (i7.empty() && i5.full())
        return true;
    // demultiplexing on i1 only
    if (i7.empty() && i5.empty() && i1.full())
        return true;

    string error_messages;
    // the following things are not allowed because they eventually dont allow assigning
    // reads unambiguously to samples
    if (i7.sparse()){
        string error_msg = "Not all samples have an i7 barcode defined. An i7 barcode needs "\
                     "to be either specified for all or none of the samples. Samples "\
                     "without a barcode: "; // + i7.samples_without_barcodes();
        error_messages.append(error_msg);
    }
    if (i5.sparse()){
	string error_msg = "Not all samples have an i5 barcode defined. An i5 barcode needs "\
                     "to be either specified for all or none of the samples. Samples "\
                     "without a barcode: "; // + i5.samples_without_barcodes();
        error_messages.append(error_msg);
    }
    if (i1.sparse()){
	string error_msg = "Not all samples have an i7 barcode defined. An i7 barcode needs "\
                     "to be either specified for all or none of the samples. Samples "\
                     "without a barcode: "; // + i1.samples_without_barcodes();
        error_messages.append(error_msg);
    }
    // an empty sample sheet does nothing and is therefore disallowed
    if (i7.empty() && i5.empty() && i1.empty()){
	string error_msg = "No index sequences have been specified in the sample sheet. "\
                     "Please specify some barcodes to run this tool.";
        error_messages.append(error_msg);
    }
    // if we did not return true earlier it is error raising time!
    fprintf(stderr, error_messages.c_str());
    return false;
}


/*
 * Generator that opens two paired fastq.gz files and yields them as utf-8 strings.

    Args:
        fq_gz_1: Path to read1 of the mate pair.
        fq_gz_2: Path to read2 of the mate pair.

    Yields:
        iter,iter: Zipped iterator for read1 and read2
*/
void Parser::get_pe_fastq(string fq_gz_1, string fq_gz_2) {
  fprintf(stderr, "TO IMPLEMENT! --> Replaced by PairedReader.cpp");
  /*
   @contextmanager
def get_pe_fastq(fq_gz_1, fq_gz_2):

    try:
        input_streams = list()
        input_streams.append(gzip.open(fq_gz_1, 'rt', encoding='utf-8'))
        input_streams.append(gzip.open(fq_gz_2, 'rt', encoding='utf-8'))
        fastq_1 = fastq_lines_to_reads(input_streams[0])
        fastq_2 = fastq_lines_to_reads(input_streams[1])
        yield zip(fastq_1, fastq_2)
    finally:
        for fq_file in input_streams:
            log.debug('Trying to close file: %s' % fq_file.name)
            try:
                fq_file.close()
            except IOError:
                log.exception("Could not close the following input file: %s"
                              % fq_file.name)

   */
}

void Parser::fastq_lines_to_reads(string fastq_lines) {
  fprintf(stderr, "TO IMPLEMENT! --> Replaced by PairedReader");
  /*
   def fastq_lines_to_reads(fastq_lines):
    """Collect fastq lines and chucks them into reads. Base on itertools
    grouper example (see https://docs.python.org/3/library/itertools.html).

    Args:
        fastq_lines (iterable <str>): An iterable of lines from a fastq file.

    Return:
        iter: All fastq lines grouped into reads (chunks of 4 lines).
    """
    # convert fastq lines into chunked iterator with a length of lines_per_read
    lines_per_read = 4
    chunked_lines = [iter(fastq_lines)] * lines_per_read
    # zip longest calls next when zipping the chunked iterator
    # this leads to all fastq_lines grouped by chunks of size lines_per_read
    return zip_longest(*chunked_lines)
   */
}

/*
 * load error correction map from misc folder.
 */
unordered_map<string,string>* Parser::get_map_from_resource(string package, string resource) {
    printf("Loading error correction map from %s", resource);
    unordered_map<string,string> *mapping = new unordered_map<string,string>();

    ifstream dataFile;
    dataFile.open(package + "/" + resource);
    while(!dataFile.eof()) {
        std::string temp;
        std::vector<string> values;
        while(std::getline(dataFile, temp, '\t') ) {
            values.push_back( string(temp.c_str()));
        }
        (*mapping)[values[0]] = values[1];
        //data.push_back(values[numColumn]);
    }
    return mapping;

}

/**
 * Reads a tsv barcodes file, converts them to byes and returns them as a dict.
    These dicts are used to map erroneous barcodes to their corrected version.
    Args:
        barcode_type (str): A path pointing to a error correction map file.
        barcode_length (int): A path pointing to a error correction map file.

    Return:
        dict (dict): correction_map : <b'erroneous_barcode', b'corrected_barcode'>
 //TODO: correction map is appended to barcode --> make this a function of the Barcode class.
 */
Barcode* Parser::load_correction_map(Barcode &barcode) {
  fprintf(stdout, "Trying to find the appropriate barcode set for {barcode.name}...\n");

    // When there are no barcodes specified, there is nothing to correct.
    if (barcode.length == 0){
        fprintf(stdout, "No barcodes have been specified for {barcode.name}.\n");
        //barcode.correction_map = NULL; //unordered_map<string,string> {"None": "None"};
        //return NULL;
    }
    unordered_set<int>* sizes = barcode.get_set_sizes();
    bool drop_none = true;
    for(auto it = sizes->begin(); it != sizes->end(); it++){
      string package_str = "./misc/barcodes/{barcode.name}";
      string file_str = "base_mapping_b{set_size}_l{barcode.length}.tsv";

      unordered_map<string,string> *corr_map = get_map_from_resource(package_str, file_str);
      //_barcode_set = set(corr_map.values());
      unordered_set<string> * _barcodes_given = barcode.get_used_codes(drop_none);
      bool is_contained = true;
		for(auto it_test = _barcodes_given->begin(); it_test != _barcodes_given->end(); it_test++){
			auto it_contained  = corr_map->find(*it_test);
			if(it_contained == corr_map->end()){
				is_contained = false;
				break;
			}
		}
		if (is_contained){
			printf("Correct set found. Used set is {set_size} barcodes with "\
                     "{barcode.length} nt length.\n");
			barcode.correction_map = *corr_map;
			return &barcode;
		}
		/*
        if _barcodes_given <= _barcode_set:
            log.info(f"Correct set found. Used set is {set_size} barcodes with "
                     f"{barcode.length} nt length.")
            barcode.correction_map = corr_map
            return barcode;
            */
    }
    if (barcode.length != 6){
    	printf("No fitting Lexogen barcode set found for {barcode.name}. No "\
                "error correction will take place for this barcode. Are you using "\
                "valid Lexogen barodes?");
    	unordered_set<string>* _bc_list = barcode.used_codes();
    	barcode.correction_map.clear();
    	for(auto it = _bc_list->begin(); it != _bc_list->end(); it++){
    		barcode.correction_map[*it] = *it;
    	}
    }
                /*
    if barcode.length != 6:
        log.warning(f"No fitting Lexogen barcode set found for {barcode.name}. No "
                    f"error correction will take place for this barcode. Are you using "
                    f"valid Lexogen barodes?")

    _bc_list = list(barcode.used_codes)
    barcode.correction_map = dict(zip(_bc_list, _bc_list))
    return barcode
   */
    return &barcode;
}

/*
    Reads the first 100 lines of paired fastq.gz files and checks if everything is
    okay with the fastq header format.

    Args:
        fq_gz_1 (str): File path of read mate1.
        fq_gz_2 (str): File path of read mate2.
        has_i7 (bool): Did the sample_sheet specify that samples have an i7 index?
        has_i5 (bool): Did the sample_sheet specify that samples have an i5 index?

    Raises:
        ValueError: When the fastq header contains less barcodes than indicated by the
            booleans.
 */
void Parser::peek_into_fastq_files(string fq_gz_1, string fq_gz_2, bool has_i7,
    bool has_i5, bool has_i1, int i7_length, int i5_length,
    int i1_start, int i1_end) {
  fprintf(stdout, "Peeking into fastq files to check for barcode formatting errors\n");

    int lines_to_check = 1000;
    int counter = 0;
    fprintf(stdout, "Checking fastq input files...\n");
    PairedReader get_pe_fastq(fq_gz_1, fq_gz_2);
    //vector<fq_read> *pe_reads = get_pe_fastq(fq_gz_1, fq_gz_2);
    std::vector<std::pair<fq_read*,fq_read*>>* pe_reads = get_pe_fastq.next_reads(lines_to_check);
        for(int i = 0; i < pe_reads->size(); i++){
          // mate_pair in pe_reads:
        	std::pair<fq_read*,fq_read*> mate_pair = pe_reads->at(i);
          check_mate_pair(mate_pair, has_i7, has_i5, has_i1, i7_length, i5_length,
                                      i1_start, i1_end);
                      counter += 1;
                      if (counter == lines_to_check)
                          break;
        }

        fprintf(stdout, "Input file formatting seems fine.\n");
}

void Parser::check_mate_pair(std::pair<fq_read*,fq_read*> mate_pair, bool has_i7, bool has_i5,
    bool has_i1, int i7_length, int i5_length,int i1_start, int i1_end) {

	fq_read* mate2 = mate_pair.second;
    try{
        check_fastq_headers(mate_pair, has_i7, has_i5, i7_length, i5_length);
        if (has_i1)
            check_mate2_length(mate2, i1_start, i1_end);
    } catch(std::exception e){
			fprintf(stderr,e.what());
			exit(EXIT_FAILURE);
    }
}

void Parser::check_mate2_length(fq_read* mate2, int i1_start, int i1_end) {
    int seq_idx = 1;
    string seq = mate2->Seq_ID; //[seq_idx]
    if (seq.length() < i1_end)
        fprintf(stderr, "Mate 2 is too short for the provided i1 barcode settings. "\
                         "According to your settings i1 starts at position {i1_start} "\
                         "and has a length of {i1_start - i1_start}. The sequence of "\
                         "mate 2 is however only {len(seq)} nt long.");

}

/**
 * Function to check if the barcodes (i7,i5) specified in the sample sheet are
    as well in the fastq header.

    Args:
        mate_pair (tuple): A tuple of mate_pairs as returned by fastq_lines_to_reads.
        has_i7 (bool): Did the sample_sheet specify that samples have an i7 index?
        has_i5 (bool): Did the sample_sheet specify that samples have an i5 index?

    Raises:
        ValueError: When the fastq header contains less barcodes than indicated by the
            booleans.
 */
void Parser::check_fastq_headers(std::pair<fq_read*,fq_read*> mate_pair, bool has_i7, bool has_i5,
    bool i7_length, int i5_length) {

    int header_idx = 0;
    fq_read* m_1 = mate_pair.first;
	fq_read* m_2 = mate_pair.second;
    //header_mate_1, header_mate_2 = m_1[header_idx], m_2[header_idx]
    string header_mate_1(m_1->Seq_ID);
    string header_mate_2(m_2->Seq_ID);

    // get the barcodes from the fastq header
    std::pair<string,string> bcs_mate1;// = header_mate_1.strip().rpartition(":")[-1].split("+")
    std::pair<string,string> bcs_mate2;// = header_mate_2.strip().rpartition(":")[-1].split("+")

    bcs_mate1 = Parser::parse_indices(header_mate_1);
    bcs_mate2 = Parser::parse_indices(header_mate_2);

    if (bcs_mate1 != bcs_mate2){
        string error_msg = "Mate1 and mate2 contain different barcode information. Please "\
                     "make sure the reads in your fastq files are paired.\n"\
                     "Mate1 header: "+header_mate_1+"\n"\
                     "Mate2 header: "+header_mate_2+"\n";
        fprintf(stderr, error_msg.c_str());
        exit(EXIT_FAILURE);
    }

    int number_bc_m1 = bcs_mate1.first == "" || bcs_mate1.second == ""  ? 1 : 2;
    int number_bc_m2 = bcs_mate2.first == "" || bcs_mate2.second == "" ? 1 : 2;

    int number_bc_present[] = {number_bc_m1, number_bc_m2};
    int expected_number = 0; //sum([has_i7, has_i5])
    if (has_i7)
    	expected_number +=1;
    if (has_i5)
    	expected_number +=2;
    // this is how a fastq header should look like
    string example_header_1 = "@NB502007:379:HM7H2BGXF:1:11101:24585:1069 1:N:0:TCAGGTAANNTT";
    string example_header_2 = "@NB502007:379:HM7H2BGXF:1:11101:24585:1069 "\
                        "1:N:0:TCAGGTAANNTT+NANGGNNCNNNN";

    // check if the header conforms to what was specified in the sample sheet
    //right_number_of_barcodes = [n <= expected_number for n in number_bc_present]
    bool right_number_of_barcodes = true;
    for(int i = 0; i < 2; i++){
    	if( number_bc_present[i] > expected_number){
    		right_number_of_barcodes = false;
    		break;
    	}
    }

    if (not right_number_of_barcodes){ //not all(right_number_of_barcodes):
        string example_header = expected_number == 2 ? example_header_2 : example_header_1;
        string error_msg = "The fastq file does not contain sufficient barcode information "\
                     "in the header.\nExpected number of barcodes: {expected_number}\n"\
                     "Observed number of barcodes: {number_bc_present}\n"\
                     "Please check your input file. Your fastq header should look "\
                     "similar to this example.\n"\
                     "Example: {example_header}\n"\
                     "Observed headers: {[header_mate_1, header_mate_2]}";
        fprintf(stderr, error_msg.c_str());
        exit(EXIT_FAILURE);
    }

    // when there are 2 barcodes in the fastq header the orientation is i7,i5
    if (has_i7 && has_i5)
        if (bcs_mate1.first.length() != i7_length || bcs_mate1.second.length() != i5_length){
		fprintf(stderr,"i7 and i5 have a different length than specified in the "\
                             "sample_sheet. "\
                             "Observed length(i7,i5): {len(bcs_mate1[0])}"\
                             ",{len(bcs_mate1[1])}\n "\
                             "Expected length(i7,i5): {i7_length},{i5_length}");
		exit(EXIT_FAILURE);
        }
    if (has_i7 && !has_i5)
        if (bcs_mate1.first.length() != i7_length){
		fprintf(stderr,"i7 has a different length than specified in the "\
                             "sample_sheet. "\
                             "Observed length(i7): {len(bcs_mate1[0])}\n"\
                             "Expected length(i7): {i7_length}\n");
		exit(EXIT_FAILURE);
        }
    if (!has_i7 && has_i5)
        if (bcs_mate1.first.length() != i5_length){
		fprintf(stderr,"i5 has a different length than specified in the "\
                             "sample_sheet. "\
                             "Observed length(i5): {len(bcs_mate1[0])}\n"\
                             "Expected length(i5): {i5_length}\n");
		exit(EXIT_FAILURE);
        }

}



Parser::~Parser() {
  // TODO Auto-generated destructor stub
}

