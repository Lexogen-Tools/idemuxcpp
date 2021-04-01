
#include "BamReader.h"

#include <iostream>
#include <fstream>
#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include "FastqReader.h"
#include <regex>
#include <string>

using namespace std;

BamReader::BamReader(string filename) {
	this->FileHandle = new BamTools::BamReader();
	bool is_open = this->FileHandle->Open(filename);
}

//const regex regex_i5("B2Z([ACGUPTN]+)");
//const regex regex_i7("BCZ([ACGUPTN]+)");
//const regex regex_bc_read1("1(\\:N\\:0\\:[ACGUPTN]*\\+?[ACGUPTN]*)$");

void append_barcodes(string &fq_header, int read_type, const string &i7, const string &i5){
    fq_header.append(" ");
    fq_header.append(to_string(read_type));
    fq_header.append(":N:0:");
    fq_header.append(i7);
    if(i7.length() > 0){
        fq_header.append("+");
    }
    fq_header.append(i5);
}

/*
   parse bam read data, which looks like this:

		 * tagdata B2ZNNNNNNNNNNNNNQ2Z#############BCZNNNNNNNNNNNNNRGZAAAHY5YM5.1QTZ#############
		VH00140:7:AAAHY5YM5:1:1101:62465:1038
		CCGGGGTCGGNGCCCCCGGCGGGCTCGGCNNCNNNNNGNNNNNNNNNNCNNNNNNNNNNNNNNNNNGNNNNNNNNNNNNN
		+
		;C;CC;CCC;#-CC;CCCCC-CC-CC;;;##C#####C##########C#################C#############

     	        <read> 	Numerical 	Read number. 1 can be single read or Read 2 of paired-end
		<is filtered> 	Y or N 	Y if the read is filtered (did not pass), N otherwise
		<control number> 	Numerical 	0 when none of the control bits are on, otherwise it is an even number. On HiSeq X systems, control specification is not performed and this number is always 0.

*/
fq_read* BamReader::get_next_read_from_bam(string &out_i7, string &out_i5, int &out_read_type){
        BamTools::BamAlignment ba;
	bool has_aln = this->FileHandle->GetNextAlignment(ba);
	if(has_aln) {
                fq_read* read = new fq_read();
		read->Seq_ID = "@" + ba.Name;

		if(ba.IsFirstMate()){
                        out_read_type = 1;
		}
		else if(ba.IsSecondMate()){
                        out_read_type = 2;
		}

                ba.GetTag("BCZ", out_i7);

                ba.GetTag("B2Z", out_i5);

		read->Sequence = ba.QueryBases;
		read->Plus_ID = "+";
		read->QualityCode = ba.Qualities;
                return read;
	}
	else {
		return NULL;
	}
}

fq_read* BamReader::next_read() {
        string i7 = "";
        string i5 = "";
        int read_type = 0;
	fq_read *read = get_next_read_from_bam(i7, i5, read_type);
        if(read != NULL) {
            append_barcodes(read->Seq_ID, read_type, i7, i5);
        }
	return read;
}

std::pair<fq_read*, fq_read*> BamReader::next_read_pair(){
    std::pair<fq_read*, fq_read*> read;

    string i7 = "";
    string i5 = "";
    int read_type = 0;
    fq_read *read1 = get_next_read_from_bam(i7, i5, read_type);
    if(read1 != NULL) {
        append_barcodes(read1->Seq_ID, read_type, i7, i5);
        // prepare read2 suffix (because i7, i5 is not contained in read 2.)
        string barcode_suffix = "";
        append_barcodes(barcode_suffix, 2, i7, i5);
        //i7 = "";
        //i5 = "";
        fq_read *read2 = get_next_read_from_bam(i7, i5, read_type);
        if(read2 != NULL){
            read2->Seq_ID.append(barcode_suffix);
            read.first = read1;
            read.second = read2;
        }
    }
    return read;
}

void BamReader::close() {
	this->FileHandle->Close();
}

BamReader::~BamReader() {
	this->close();
	delete this->FileHandle;
}

