#ifndef BARCODE_H_
#define BARCODE_H_

#include <vector>
#include <unordered_set>
#include <unordered_map>
using namespace std;

class Barcode {
public:
  Barcode(string barcode_type,  unordered_map<string,std::vector<string>> &ix_barcodes, bool reverse_complement = false);
  vector<int> Lengths;

  // attribute getters
  inline
  vector<string>* samples_without_barcodes(){
	  vector<string>* keys = new vector<string>();
	  for(auto it = _sample_map.begin(); it != _sample_map.end(); it++){
		  if (it->first.compare("") == 0){
			  keys->push_back(it->first);
		  }
	  }
	  return keys;
  };
  inline
  unordered_set<string>*used_codes(){
	  unordered_set<string>* keys = new unordered_set<string>();
	  for(auto it = _sample_map.begin(); it != _sample_map.end(); it++){
		  keys->insert(it->first);
	  }
	  return keys;
  };
  inline
    bool sparse(){
	  unordered_set<string>* uc = this->used_codes();
	  bool is_sparse = false;
	  for(auto it = uc->begin(); it != uc->end(); it++){
		  if (string("").compare(*it) == 0){
			  is_sparse = true;
			  break;
		  }
	  }
	  delete uc;
	  return is_sparse;
    };
  inline
    bool empty(){
	  unordered_set<string>* uc = this->used_codes();
	  //None in self.used_codes and len(self.used_codes) == 1
	  bool is_empty = (this->sparse() && uc->size() == 1) || uc->size() == 0;
	  delete uc;
	  return is_empty;
    };
  inline
    bool full(){
	  unordered_set<string>* uc = this->used_codes();
	  bool is_full = true;
	  for(auto it = uc->begin(); it != uc->end(); it++){
		  if (string("").compare(*it) == 0){
			  is_full = false;
			  break;
		  }
	  }
	  delete uc;
	  return is_full;
    };

  inline
  unordered_set<string>* get_used_codes(bool drop_none=false){
	  unordered_set<string>* _return_val = this->used_codes();
      if (drop_none){
    	  _return_val->erase("");
    	  /*
    	  vector<string> to_erase;
    	  for(auto it = _return_val->begin(); it != _return_val->end(); it++){
    		  if(string("").compare(*it) == 0){
    			  to_erase.push_back(*it);
    		  }
    		  else{
    			  i++;
    		  }
    	  }
    	  */
         // _return_val.discard(None);
      }
      return _return_val;
  }

  void post_init(bool reverse_complement);
  void check_length();
  void load_correction_map(string relative_exepath, string correction_maps_path = "");
  void create_one_to_one_map();

  //barcode name is i7,i5,i1 (name of the resources directory).
  string Name;
  unordered_map<string,string> Correction_map;
  unordered_map<string,string> _sample_map;
  unordered_set<int> allowed_lengths;
  unordered_set<int> lengths_96_barcodes;
  unordered_set<int> lengths_384_barcodes;
  unordered_set<int>* get_set_sizes();
  virtual ~Barcode();
};

#endif /* BARCODE_H_ */


