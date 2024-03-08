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


template<typename ... Args> string
string_format(string  format,
              Args... args)
{
  char  buffer[100000];
  int   cx = snprintf(buffer, 100000, format.c_str(), args ...);

  if (cx < 0)
    throw(std::runtime_error("Error: could not format the input pattern! " + format));

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
static bool
folder_exists(std::string foldername)
{
  boost::filesystem::path p(foldername);

  return boost::filesystem::exists(p);
}


/**
 * Recursive, portable wrapper for mkdir.
 * @param[in] path the full path of the directory to create.
 * @return zero on success, otherwise -1.
 */
static int
mkdir(string path)
{
  boost::filesystem::path tmp_path(path);

  return boost::filesystem::create_directory(path) ? 0 : -1;
}


static int
rmdir(string path)
{
  return boost::filesystem::remove_all(path) ? 0 : -1;
}


/*
 * Returns the full path to the currently running executable,
 * or an empty string in case of failure.
 */
static std::string
getExecutablePath()
{
#if (BOOST_OS_WINDOWS)
  char    *exePath;
  if (_get_pgmptr(&exePath) != 0)
    exePath = "";

#elif (BOOST_OS_SOLARIS)
  char    exePath[PATH_MAX];
  if (realpath(getexecname(), exePath) == NULL)
    exePath[0] = '\0';

#elif (BOOST_OS_LINUX)
  char    exePath[PATH_MAX];
  ssize_t len = ::readlink("/proc/self/exe", exePath, sizeof(exePath));
  if (len == -1 || len == sizeof(exePath))
    len = 0;

  exePath[len] = '\0';
#elif (BOOST_OS_MACOS)
  char      exePath[PATH_MAX];
  uint32_t  len = sizeof(exePath);
  if (_NSGetExecutablePath(exePath, &len) != 0) {
    exePath[0] = '\0';     // buffer too small (!)
  } else {
    // resolve symlinks, ., .. if possible
    char *canonicalPath = realpath(exePath, NULL);
    if (canonicalPath != NULL) {
      strncpy(exePath, canonicalPath, len);
      free(canonicalPath);
    }
  }

#elif (BOOST_OS_BSD_FREE)
  char    exePath[2048];
  int     mib[4];
  mib[0]  = CTL_KERN;
  mib[1]  = KERN_PROC;
  mib[2]  = KERN_PROC_PATHNAME;
  mib[3]  = -1;
  size_t  len = sizeof(exePath);
  if (sysctl(mib, 4, exePath, &len, NULL, 0) != 0)
    exePath[0] = '\0';

#endif
  return std::string(exePath);
}
};

#endif /* SRC_HELPER_H_ */
