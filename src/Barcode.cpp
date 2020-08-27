/*
 * Barcode.cpp
 *
 *  Created on: Aug 15, 2020
 *      Author: quin
 */

#include <vector>
#include <string>
#include <unordered_map>
#include "Barcode.h"

using namespace std;

Barcode::Barcode(string barcode_type,  unordered_map<string,std::vector<string>> &ix_barcodes, bool reverse_complement) : length(0),
    lengths_96_barcodes({8, 10, 12}),
    lengths_384_barcodes({10, 12}), Name("") {
  // TODO Auto-generated constructor stub
	for(auto it = ix_barcodes.begin(); it != ix_barcodes.end(); it++)
		this->_sample_map[it->first] = it->second[0];
	post_init(reverse_complement);
}

unordered_set<int>* Barcode::get_set_sizes(){
  unordered_set<int>* set_sizes = new unordered_set<int>();
  auto it = allowed_lengths.find(length);
  if (it != allowed_lengths.end()){
      auto it_96 = lengths_96_barcodes.find(length);
      if(it_96 != lengths_96_barcodes.end())
          set_sizes->insert(96);
      auto it_384 = lengths_384_barcodes.find(length);
      if(it_384 != lengths_384_barcodes.end())
          set_sizes->insert(384);
  }
  return set_sizes;
}


void Barcode::post_init(bool reverse_complement){
    // if the barcode is reverse complement add _rc to the name
    // we need this later to load the correction maps
    if (reverse_complement)
        this->Name += "_rc";
    this->check_length();
}

void Barcode::check_length(){
// barcodes for each type need to be of the same length. This is as the
// sequences of 8 nt long barcodes are contained within 10 nt long barcodes,
// and the ones of 10 in 12. if we dont enforce the same length per barcode there
// might be a possibility we cant tell with certainty to which sample a
// barcoded read belongs.
	for(auto it = this->_sample_map.begin(); it != this->_sample_map.end(); it++){
		string barcode = it->first;
		string sample_name = it->second;
		if (barcode.compare("") != 0){
			auto itl = this->allowed_lengths.find(barcode.length());
			if(itl != this->allowed_lengths.end()){
				this->observed_lengths.insert(barcode.length());
				this->length = barcode.length();
				if(this->observed_lengths.size() > 1){
					fprintf(stderr, "{self.name} barcodes with a different length "\
                                 "have been observed for {sample_name}. Barcodes"\
                                 " need to have the same length for all "\
                                 "samples.\nObserved barcode length:"\
                                 " {len(barcode)} \nPrevious observed length:"\
                                 " {self._observed_lengths.pop()}");
					exit(EXIT_FAILURE);
				}
			}
			else{
				fprintf(stderr,"{self.name} barcodes are {len(barcode)} nt long."\
                             "{self.name} barcodes are only allowed to be"\
                             " {self._allowed_lengths} nt long.");
				exit(EXIT_FAILURE);
			}

		}
	}
}
/*


@dataclass
class Barcode:
    name: str
    _sample_map: defaultdict = field(default_factory=lambda: defaultdict(list))
    reverse_complement: bool = False
    _observed_lengths: set = field(default_factory=set)
    correction_map: dict = None
    length: int = None
    _lengths_96_barcodes = {8, 10, 12}
    _lengths_384_barcodes = {10, 12}
    _allowed_lengths = {6, 8, 10, 12}

    def __post_init__(self):
        # if the barcode is reverse complement add _rc to the name
        # we need this later to load the correction maps
        if self.reverse_complement:
            self.name = f"{self.name}_rc"
        self._check_length()

    def _check_length(self):
        # barcodes for each type need to be of the same length. This is as the
        # sequences of 8 nt long barcodes are contained within 10 nt long barcodes,
        # and the ones of 10 in 12. if we dont enforce the same length per barcode there
        # might be a possibility we cant tell with certainty to which sample a
        # barcoded read belongs.
        for barcode, sample_name in self._sample_map.items():
            if barcode is not None:
                if len(barcode) in self._allowed_lengths:
                    self._observed_lengths.add(len(barcode))
                    self.length = len(barcode)
                    if len(self._observed_lengths) > 1:
                        raise ValueError(f"{self.name} barcodes with a different length "
                                         f"have been observed for {sample_name}. Barcodes"
                                         f" need to have the same length for all "
                                         f"samples.\nObserved barcode length:"
                                         f" {len(barcode)} \nPrevious observed length:"
                                         f" {self._observed_lengths.pop()}")
                else:
                    raise ValueError(f"{self.name} barcodes are {len(barcode)} nt long."
                                     f"{self.name} barcodes are only allowed to be"
                                     f" {self._allowed_lengths} nt long.")

    def get_set_sizes(self):
        set_sizes = []
        if self.length in self._allowed_lengths:
            if self.length in self._lengths_96_barcodes:
                set_sizes.append(96)
            if self.length in self._lengths_384_barcodes:
                set_sizes.append(384)
        return set_sizes

    def get_used_codes(self, drop_none=False):
        _return_val = self.used_codes.copy()
        if drop_none:
            _return_val.discard(None)
        return _return_val

    # attribute getters
    @property
    def samples_without_barcodes(self):
        return self._sample_map.get(None)

    @property
    def used_codes(self):
        return set(self._sample_map.keys())

    @property
    def sparse(self):
        return None in self.used_codes and len(self.used_codes) > 1

    @property
    def empty(self):
        return None in self.used_codes and len(self.used_codes) == 1

    @property
    def full(self):
        return None not in self.used_codes
 */

Barcode::~Barcode() {
  // TODO Auto-generated destructor stub
}

