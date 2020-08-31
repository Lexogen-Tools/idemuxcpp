#include <string>
#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include "ZipFastqReader.h"

ZipFastqReader::ZipFastqReader(string filename) : BUFLEN(8024), buf(new char[BUFLEN]),offset(buf),cur(NULL),end(NULL),read_line_index(0) {
	GZ_filehandle = gzopen(filename.c_str(), "rb");
	if(GZ_filehandle == NULL)
		std::cerr << "Error: could not open gzip file: " << filename << std::endl;
}

void error(const char* const msg)
{
    std::cerr << msg << "\n";
    exit(255);
}

fq_read* ZipFastqReader::next_read(){
    //fq_read* read = read_in;
    //if(!read_in)
    //	read = new fq_read();
	fq_read* read = new fq_read();
    int err, len;
    read_line_index = 0;
    // finish prev. readbuffer.
    if (cur < end)
    	goto label_printbuffer;
    // start new read.
    for (;;) {
        //int err, len = sizeof(buf)-(offset-buf);
        len = BUFLEN-(offset-buf);

        if (len == 0) error("Buffer to small for input line lengths");

        len = gzread(GZ_filehandle, offset, len);

        if (len == 0){
        	break;
        }
        if (len <  0) error(gzerror(GZ_filehandle, &err));

        cur = buf;
        end = offset+len;
        label_printbuffer:
        for (char* eol; (cur<end) && (eol = std::find(cur, end, '\n')) < end; )
        {
        	switch(read_line_index){
        	case 0:
        		//read->Seq_ID = std::string(cur, eol);
        		read->Seq_ID = string(cur,eol-cur);
        		break;
        	case 1:
				//read->Sequence = std::string(cur, eol);
        		read->Sequence = string(cur,eol-cur);
				break;
        	case 2:
				//read->Plus_ID = std::string(cur, eol);
        		read->Plus_ID = string(cur,eol-cur);
				break;
        	case 3:
				//read->QualityCode = std::string(cur, eol);
        		read->QualityCode = string(cur,eol-cur);
				break;

        	}
        	read_line_index +=1;
        	if(read_line_index == 4){
        		cur = eol + 1;
        		return read;
        	}
        	cur = eol + 1;
        }

        // any trailing data in [eol, end) now is a partial line
        offset = std::copy(cur, end, buf);
        //offset = strncpy(offset,cur,)
    }

    // BIG CATCH: don't forget about trailing data without eol :)
    //std::cout << std::string(buf, offset);
    if ((cur<end) && (!std::find(cur, end, '\n'))){
		read_line_index +=1;
		switch(read_line_index){
		case 0:
			//read->Seq_ID = std::string(cur, eol);
			read->Seq_ID = string(buf,buf-offset);
			break;
		case 1:
			//read->Sequence = std::string(cur, eol);
			read->Sequence = string(buf,buf-offset);
			break;
		case 2:
			//read->Plus_ID = std::string(cur, eol);
			read->Plus_ID = string(buf,buf-offset);
			break;
		case 3:
			//read->QualityCode = std::string(cur, eol);
			read->QualityCode = string(buf,buf-offset);
			break;

		}
    }
    else{
    	if (len == 0){
    		delete read;
    		//read->Seq_ID = "";
    		//return read;
    		return NULL;
    	}
    }

	return read;
}

void ZipFastqReader::close(){
	if (gzclose(GZ_filehandle) != Z_OK) error("failed gzclose");
}


ZipFastqReader::~ZipFastqReader() {
	close();
	delete[] buf;
}

