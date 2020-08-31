
#ifndef SRC_HELPER_H_
#define SRC_HELPER_H_
#include <cstdarg>
#include <string>
#include <stdexcept>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

using namespace std;


template<typename... Args> string string_format(string format, Args... args){
char buffer [100000];
int cx = snprintf(buffer, 100000, format.c_str(), args...);
if(cx < 0)
	throw(std::runtime_error("Error: could not format the input pattern! "+ format));
return string(buffer);
}

class utils
{
public:
    /**
     * Checks if a folder exists
     * @param foldername path to the folder to check.
     * @return true if the folder exists, false otherwise.
     */
    static bool folder_exists(std::string foldername)
    {
    	boost::filesystem::path p(foldername);
    	return boost::filesystem::exists(p);
    }

    /**
     * Recursive, portable wrapper for mkdir.
     * @param[in] path the full path of the directory to create.
     * @return zero on success, otherwise -1.
     */
    static int mkdir(string path)
    {
        boost::filesystem::path tmp_path(path);
        return boost::filesystem::create_directory(path) ? 0 : -1;
    }

    static int rmdir(string path){
    	return boost::filesystem::remove_all(path) ? 0 : -1;
    }
};

#endif /* SRC_HELPER_H_ */
