/** @file idemuxCPP_cmdl.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.6
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef IDEMUXCPP_CMDL_H
#define IDEMUXCPP_CMDL_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef IDEMUXCPP_CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define IDEMUXCPP_CMDLINE_PARSER_PACKAGE "idemuxCPP"
#endif

#ifndef IDEMUXCPP_CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#define IDEMUXCPP_CMDLINE_PARSER_PACKAGE_NAME "idemuxCPP"
#endif

#ifndef IDEMUXCPP_CMDLINE_PARSER_VERSION
/** @brief the program version */
#define IDEMUXCPP_CMDLINE_PARSER_VERSION VERSION
#endif

/** @brief Where the command line options are stored */
struct idemuxCPP_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  char * r1_arg;	/**< @brief Fastq.gz read file 1.
 (default='').  */
  char * r1_orig;	/**< @brief Fastq.gz read file 1.
 original value given at command line.  */
  const char *r1_help; /**< @brief Fastq.gz read file 1.
 help description.  */
  char * r2_arg;	/**< @brief Fastq.gz read file 2.
 (default='').  */
  char * r2_orig;	/**< @brief Fastq.gz read file 2.
 original value given at command line.  */
  const char *r2_help; /**< @brief Fastq.gz read file 2.
 help description.  */
  char * out_arg;	/**< @brief Where to write the output files.
 (default='./').  */
  char * out_orig;	/**< @brief Where to write the output files.
 original value given at command line.  */
  const char *out_help; /**< @brief Where to write the output files.
 help description.  */
  char * sample_sheet_arg;	/**< @brief Outputs a csv file containing sample names, i7, i5 and i1 barcodes.
 (default='sample-sheet.csv').  */
  char * sample_sheet_orig;	/**< @brief Outputs a csv file containing sample names, i7, i5 and i1 barcodes.
 original value given at command line.  */
  const char *sample_sheet_help; /**< @brief Outputs a csv file containing sample names, i7, i5 and i1 barcodes.
 help description.  */
  char * barcode_corrections_arg;	/**< @brief Outputs a csv file that contains the number of corrected barcodes.  */
  char * barcode_corrections_orig;	/**< @brief Outputs a csv file that contains the number of corrected barcodes original value given at command line.  */
  const char *barcode_corrections_help; /**< @brief Outputs a csv file that contains the number of corrected barcodes help description.  */
  int i5_rc_flag;	/**< @brief Set this flag if the i5 barcode has been sequenced as reverse complement and the barcodes you provided should be reverse complemented.
 (default=off).  */
  const char *i5_rc_help; /**< @brief Set this flag if the i5 barcode has been sequenced as reverse complement and the barcodes you provided should be reverse complemented.
 help description.  */
  int i1_start_arg;	/**< @brief Start position of the i1 index (1-based) on read 2.
 (default='11').  */
  char * i1_start_orig;	/**< @brief Start position of the i1 index (1-based) on read 2.
 original value given at command line.  */
  const char *i1_start_help; /**< @brief Start position of the i1 index (1-based) on read 2.
 help description.  */
  int queue_size_arg;	/**< @brief Queue size for reads that will be processed in one block.
 (default='4000000').  */
  char * queue_size_orig;	/**< @brief Queue size for reads that will be processed in one block.
 original value given at command line.  */
  const char *queue_size_help; /**< @brief Queue size for reads that will be processed in one block.
 help description.  */
  int reading_threads_arg;	/**< @brief Number of threads used for reading gz files. Either 1 or 2 (one thread per input file is used).
 (default='2').  */
  char * reading_threads_orig;	/**< @brief Number of threads used for reading gz files. Either 1 or 2 (one thread per input file is used).
 original value given at command line.  */
  const char *reading_threads_help; /**< @brief Number of threads used for reading gz files. Either 1 or 2 (one thread per input file is used).
 help description.  */
  int writing_threads_arg;	/**< @brief Number of threads used for writing gz files. Default is the number of processor cores.
.  */
  char * writing_threads_orig;	/**< @brief Number of threads used for writing gz files. Default is the number of processor cores.
 original value given at command line.  */
  const char *writing_threads_help; /**< @brief Number of threads used for writing gz files. Default is the number of processor cores.
 help description.  */
  int processing_threads_arg;	/**< @brief Number of threads used for processing the error correction. Default is the number of processor cores.
.  */
  char * processing_threads_orig;	/**< @brief Number of threads used for processing the error correction. Default is the number of processor cores.
 original value given at command line.  */
  const char *processing_threads_help; /**< @brief Number of threads used for processing the error correction. Default is the number of processor cores.
 help description.  */
  int demux_only_flag;	/**< @brief Do a one on one mapping for the barcodes specified in the sample sheet. No error correction will be done. Barcodes that do not match are written to the undetermined reads file. (default=off).  */
  const char *demux_only_help; /**< @brief Do a one on one mapping for the barcodes specified in the sample sheet. No error correction will be done. Barcodes that do not match are written to the undetermined reads file. help description.  */
  int verbose_flag;	/**< @brief Verbose.
 (default=off).  */
  const char *verbose_help; /**< @brief Verbose.
 help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int r1_given ;	/**< @brief Whether r1 was given.  */
  unsigned int r2_given ;	/**< @brief Whether r2 was given.  */
  unsigned int out_given ;	/**< @brief Whether out was given.  */
  unsigned int sample_sheet_given ;	/**< @brief Whether sample-sheet was given.  */
  unsigned int barcode_corrections_given ;	/**< @brief Whether barcode-corrections was given.  */
  unsigned int i5_rc_given ;	/**< @brief Whether i5-rc was given.  */
  unsigned int i1_start_given ;	/**< @brief Whether i1-start was given.  */
  unsigned int queue_size_given ;	/**< @brief Whether queue-size was given.  */
  unsigned int reading_threads_given ;	/**< @brief Whether reading-threads was given.  */
  unsigned int writing_threads_given ;	/**< @brief Whether writing-threads was given.  */
  unsigned int processing_threads_given ;	/**< @brief Whether processing-threads was given.  */
  unsigned int demux_only_given ;	/**< @brief Whether demux-only was given.  */
  unsigned int verbose_given ;	/**< @brief Whether verbose was given.  */

} ;

/** @brief The additional parameters to pass to parser functions */
struct idemuxCPP_cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure idemuxCPP_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure idemuxCPP_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *idemuxCPP_args_info_purpose;
/** @brief the usage string of the program */
extern const char *idemuxCPP_args_info_usage;
/** @brief the description string of the program */
extern const char *idemuxCPP_args_info_description;
/** @brief all the lines making the help output */
extern const char *idemuxCPP_args_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int idemuxCPP_cmdline_parser (int argc, char **argv,
  struct idemuxCPP_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use idemuxCPP_cmdline_parser_ext() instead
 */
int idemuxCPP_cmdline_parser2 (int argc, char **argv,
  struct idemuxCPP_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int idemuxCPP_cmdline_parser_ext (int argc, char **argv,
  struct idemuxCPP_args_info *args_info,
  struct idemuxCPP_cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int idemuxCPP_cmdline_parser_dump(FILE *outfile,
  struct idemuxCPP_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int idemuxCPP_cmdline_parser_file_save(const char *filename,
  struct idemuxCPP_args_info *args_info);

/**
 * Print the help
 */
void idemuxCPP_cmdline_parser_print_help(void);
/**
 * Print the version
 */
void idemuxCPP_cmdline_parser_print_version(void);

/**
 * Initializes all the fields a idemuxCPP_cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void idemuxCPP_cmdline_parser_params_init(struct idemuxCPP_cmdline_parser_params *params);

/**
 * Allocates dynamically a idemuxCPP_cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized idemuxCPP_cmdline_parser_params structure
 */
struct idemuxCPP_cmdline_parser_params *idemuxCPP_cmdline_parser_params_create(void);

/**
 * Initializes the passed idemuxCPP_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void idemuxCPP_cmdline_parser_init (struct idemuxCPP_args_info *args_info);
/**
 * Deallocates the string fields of the idemuxCPP_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void idemuxCPP_cmdline_parser_free (struct idemuxCPP_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int idemuxCPP_cmdline_parser_required (struct idemuxCPP_args_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* IDEMUXCPP_CMDL_H */
