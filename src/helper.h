
#ifndef SRC_HELPER_H_
#define SRC_HELPER_H_
#include <cstdarg>
#include <string>
#include <stdexcept>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

//for determining the executable path (unit tests don't have argv[0])
#include <boost/predef/os.h>
#if (BOOST_OS_WINDOWS)
  #include <stdlib.h>
#elif (BOOST_OS_SOLARIS)
  #include <stdlib.h>
  #include <limits.h>
#elif (BOOST_OS_LINUX)
  #include <unistd.h>
  #include <limits.h>
#elif (BOOST_OS_MACOS)
  #include <mach-o/dyld.h>
#elif (BOOST_OS_BSD_FREE)
  #include <sys/types.h>
  #include <sys/sysctl.h>
#endif



// _setmaxstdio(2048); is in stdio.h on windows.
#ifdef _WIN32
  #include <stdio.h>
  #include <stdlib.h>
#else
  #include <sys/time.h>
  #include <sys/resource.h>
#endif


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

/*
 * Returns the full path to the currently running executable,
 * or an empty string in case of failure.
 */
static std::string getExecutablePath() {
#if (BOOST_OS_WINDOWS)
    char *exePath;
    if (_get_pgmptr(&exePath) != 0)
        exePath = "";
#elif (BOOST_OS_SOLARIS)
    char exePath[PATH_MAX];
    if (realpath(getexecname(), exePath) == NULL)
        exePath[0] = '\0';
#elif (BOOST_OS_LINUX)
    char exePath[PATH_MAX];
    ssize_t len = ::readlink("/proc/self/exe", exePath, sizeof(exePath));
    if (len == -1 || len == sizeof(exePath))
        len = 0;
    exePath[len] = '\0';
#elif (BOOST_OS_MACOS)
    char exePath[PATH_MAX];
    uint32_t len = sizeof(exePath);
    if (_NSGetExecutablePath(exePath, &len) != 0) {
        exePath[0] = '\0'; // buffer too small (!)
    } else {
        // resolve symlinks, ., .. if possible
        char *canonicalPath = realpath(exePath, NULL);
        if (canonicalPath != NULL) {
            strncpy(exePath,canonicalPath,len);
            free(canonicalPath);
        }
    }
#elif (BOOST_OS_BSD_FREE)
    char exePath[2048];
    int mib[4];  mib[0] = CTL_KERN;  mib[1] = KERN_PROC;  mib[2] = KERN_PROC_PATHNAME;  mib[3] = -1;
    size_t len = sizeof(exePath);
    if (sysctl(mib, 4, exePath, &len, NULL, 0) != 0)
        exePath[0] = '\0';
#endif
    return std::string(exePath);
}

/*
 * return true if the given limit is smaller, than the sytem limit or if it has been increased.
 * return false if the limit could not be requested or set.
 */
static
bool set_maximal_file_limit(size_t new_max_limit){
	//check if the number of barcodes (corresponds to the number of files) is higher than the maximum limit of open file handles, allowed by the system.
#ifndef _WIN32
	rlimit old_limit, new_limit;
	int received_limit = getrlimit(RLIMIT_NOFILE, &old_limit);
	if(received_limit == 0)
	{
		printf("Open files limit soft %ld, hard %ld, required %ld\n", old_limit.rlim_cur, old_limit.rlim_max, new_max_limit);
		if((old_limit.rlim_cur < new_max_limit) || (old_limit.rlim_max < new_max_limit)){
			size_t real_new_limit = min(new_max_limit,(size_t)RLIM_INFINITY);
			if (old_limit.rlim_cur < new_max_limit)
			    new_limit.rlim_cur = real_new_limit;
			if (old_limit.rlim_max < new_max_limit)
			    new_limit.rlim_max = real_new_limit;
			int set_limit = setrlimit(RLIMIT_NOFILE, &new_limit);
			if(set_limit == 0){
				fprintf(stdout, "Maximum open file limit has been increased to %ld\n", real_new_limit);
				return true;
			}
			else{
				fprintf(stderr, "Error: the maximum open file limit for the filesystem could not be set!\n");
				return false;
			}
		}
		else{
			//increasing was not necessary
			return true;
		}
	}
	else{
		fprintf(stderr, "Error: the maximum open file limit for the filesystem could not be retrieved!\n");
		return false;
	}
#else
	int current_limit = _getmaxstdio();
	if(current_limit < new_max_limit){
		printf("Open files limit %d, new %zu\n", current_limit, new_max_limit);
		int set_limit = _setmaxstdio(new_max_limit);
		if(!set_limit){
			fprintf(stderr, "Error: the maximum file limit for the filesystem could not be set!\n");
			return false;
		}
		fprintf(stdout, "Maximum open file limit has been increased to %zu\n", new_max_limit);
		return true;
	}
	else{
		//increasing was not necessary
		return true;
	}

#endif
	return false;
}
};


#endif /* SRC_HELPER_H_ */
