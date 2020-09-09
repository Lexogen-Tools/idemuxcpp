#define BOOST_TEST_MODULE boost_test_barcode
#include <boost/test/included/unit_test.hpp>
#include <unordered_map>
#include <string>
#include "../src/Barcode.h"
#include <stdexcept>

using namespace std;


unordered_map<string, std::vector<string>>* invalid_lengths(){
	unordered_map<string, std::vector<string>>* vals = new unordered_map<string, std::vector<string>>();
	vals->insert({string(7,'A'), {"test0"}});
	vals->insert({string(7,'C'), {"test1"}});
    return vals;
}

unordered_map<string, std::vector<string>>* none_values(){
	unordered_map<string, std::vector<string>>* vals = new unordered_map<string, std::vector<string>>();
	vals->insert({"", {"test0"}});
    return vals;
}

unordered_map<string, std::vector<string>>* mixed_none_values(){
	unordered_map<string, std::vector<string>>* vals = new unordered_map<string, std::vector<string>>();
	vals->insert({"", {"test0"}});
	vals->insert({string(10,'C'), {"test1"}});
    return vals;
}

unordered_map<string, std::vector<string>>* different_lengths(){
	unordered_map<string, std::vector<string>>* vals = new unordered_map<string,  std::vector<string>>();
	vals->insert({string(10,'A'), {"test0"}});
	vals->insert({string(12,'C'), {"test1"}});
    return vals;
}

unordered_map<string, std::vector<string>>* correct_length(){
	unordered_map<string, std::vector<string>>* vals = new unordered_map<string, std::vector<string>>();
	vals->insert({string(12,'A'), {"test0"}});
    return vals;
}


BOOST_AUTO_TEST_CASE( test_invalid_lengths )
{
	unordered_map<string,  std::vector<string>>* l = invalid_lengths();
	BOOST_CHECK_THROW( Barcode b("invalid length", *l,false), std::runtime_error);
	delete l;
}

BOOST_AUTO_TEST_CASE( test_different_lengths )
{
	unordered_map<string,  std::vector<string>>* l = different_lengths();
	BOOST_CHECK_THROW( Barcode b("different lengths", *l,false), std::runtime_error);
	delete l;
}

BOOST_AUTO_TEST_CASE( test_barcode_rc )
{
	unordered_map<string,  std::vector<string>>* l = correct_length();
	Barcode bc("i5", *l, true);
	BOOST_CHECK(string("i5_rc").compare(bc.Name) == 0);
	delete l;
}

BOOST_AUTO_TEST_CASE( test_is_full )
{
	unordered_map<string,  std::vector<string>>* l = correct_length();
	Barcode bc("i5", *l);
	BOOST_CHECK(bc.full() == true);
	delete l;
}

BOOST_AUTO_TEST_CASE( test_is_sparse )
{
	unordered_map<string,  std::vector<string>>* l = mixed_none_values();
	Barcode bc("i5", *l);
	BOOST_CHECK(bc.sparse() == true);
	delete l;
}

BOOST_AUTO_TEST_CASE( test_is_empty )
{
	unordered_map<string,  std::vector<string>>* l = none_values();
	Barcode bc("i5", *l);
	BOOST_CHECK(bc.empty() == true);
	BOOST_CHECK(bc.length == 0);
	delete l;
}

BOOST_AUTO_TEST_CASE( test_length_getter )
{
	unordered_map<string,  std::vector<string>>* l = correct_length();
	Barcode bc("i5", *l);
	BOOST_CHECK(bc.length == 12);
	delete l;
}



