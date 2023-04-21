
#include <thread>
#include <math.h>
#include "config.h"
#include "PairedReader.h"
#include "FastqReader.h"
#include "BoostZipReader.h"
#if HAVE_LIBBAMTOOLS
#include "BamReader.h"
#endif

PairedReader::PairedReader(string fastqgz_1, string fastqgz_2, double queue_buffer_gb) : Queue_buffer_bytes(size_t(queue_buffer_gb*pow(1024,3))){

	if (fastqgz_1.length() > 3){
		string suffix = fastqgz_1.substr(fastqgz_1.length()-3);
		if (suffix.compare(".gz") == 0){
			// use zip reader
			Reader1 = new BoostZipReader(fastqgz_1);
		}
		else{
			// use fastq reader
			Reader1 = new FastqReader(fastqgz_1);
		}
	}
	if (fastqgz_2.length() > 3){
		string suffix = fastqgz_2.substr(fastqgz_2.length()-3);
		if (suffix.compare(".gz") == 0){
			// use zip reader
			Reader2 = new BoostZipReader(fastqgz_2);
		}
		else{
			// use fastq reader
			Reader2 = new FastqReader(fastqgz_2);
		}
	}
	IsBam = false;
}

PairedReader::PairedReader(string bam_file_alternating, double queue_buffer_gb) : Queue_buffer_bytes(size_t(queue_buffer_gb*pow(1024,3))) {

	if (bam_file_alternating.length() > 0){
#if HAVE_LIBBAMTOOLS
		Reader1 = new BamReader(bam_file_alternating);
#else
		throw runtime_error("Error: bam files are not supported with this compilation!");
#endif
	}
	IsBam = true;
	Reader2 = NULL;
}

void PairedReader::fill_reads_interleaved(IFastqReader* reader, size_t max_size, std::vector<std::pair<fq_read*, fq_read*>> *pairs) {
	size_t bytes_read = 0;
	while((Queue_buffer_bytes > 0 && (bytes_read < Queue_buffer_bytes)) ||
			(max_size > 0 && (pairs->size() < max_size))){
		std::pair<fq_read*, fq_read*> pair = reader->next_read_pair();
		//pair.first = reader->next_read();
		//pair.second = reader->next_read();
		if(pair.first != NULL && pair.second != NULL){
			pairs->push_back(pair);
			if(Queue_buffer_bytes > 0){
				bytes_read += pair.first->get_size();
				bytes_read += pair.second->get_size();
			}
		}
		else{
			break;
		}
	}
}

std::vector<std::pair<fq_read*, fq_read*>>* PairedReader::next_reads(
		size_t max_size) {
	std::vector<std::pair<fq_read*, fq_read*>> *pairs = new std::vector<std::pair<fq_read*, fq_read*>>();
	//fq_read* pairs = new fq_read[max_size*2+1];
	if(IsBam){
		fill_reads_interleaved(Reader1, max_size, pairs);
	}
	else{
		fq_read *read1, *read2;
		size_t bytes_read = 0;
		while((Queue_buffer_bytes > 0 && (bytes_read < Queue_buffer_bytes)) ||
				(max_size > 0 && (pairs->size() < max_size))){
			read1 = Reader1->next_read();
			read2 = Reader2->next_read();
			if (read1 != NULL && read2 != NULL) {
				pairs->push_back(std::pair<fq_read*, fq_read*>(read1, read2));
				if(Queue_buffer_bytes > 0){
					bytes_read += read1->get_size();
					bytes_read += read2->get_size();
				}
			}
			else{
				delete read1;
				delete read2;
				break;
			}
		}
	}
	return pairs;
}

void fill_reads(IFastqReader* reader, size_t max_size, fq_read ***reads) {
	for (size_t i = 0; i < max_size; i++) {
		fq_read *tmp = reader->next_read();
		if (tmp != NULL){
			(*reads)[i] = tmp;
		}
		else{
			(*reads)[i] = NULL;
			delete tmp;
			break;
		}
	}
}

std::vector<std::pair<fq_read*, fq_read*>>* PairedReader::next_reads2(
		size_t max_size, int reading_threads) {
	std::vector<std::pair<fq_read*, fq_read*>> *pairs = new std::vector<std::pair<fq_read*, fq_read*>>();

	if(IsBam){
		fill_reads_interleaved(Reader1, max_size, pairs);
	}
	else{
		if (reading_threads == 2){
			fq_read** reads1 = new fq_read*[max_size];
			fq_read** reads2 = new fq_read*[max_size];

			std::thread to1(fill_reads,  Reader1, max_size, &reads1);
			std::thread to2(fill_reads,  Reader2, max_size, &reads2);
			to1.join();
			to2.join();

			for(size_t i = 0; i < max_size; i++){
				if(reads1[i] != NULL && reads2[i] != NULL){
					pairs->push_back(std::pair<fq_read*, fq_read*>(reads1[i], reads2[i]));
				}
				else{
					break;
				}
			}
			delete[] reads1;
			delete[] reads2;
		}
		else{
			// use 1 thread
			pairs = this->next_reads(max_size);
		}
	}
	return pairs;
}

PairedReader::~PairedReader() {
	delete Reader1;
	delete Reader2;
}

