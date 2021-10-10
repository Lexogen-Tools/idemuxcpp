AUTOMAKE_OPTIONS = subdir-objects
TESTS = test_barcode test_demuxer test_parser
check_PROGRAMS = test_barcode test_demuxer test_parser
test_barcode_CPPFLAGS = -I../src/ ${BOOST_CPPFLAGS}
test_barcode_LDADD = ${BOOST_LDFLAGS} ${LIBBOOST_UNIT_TEST_FRAMEWORK} $(LIBBOOST_FILESYSTEM_LIB) $(LIBBOOST_SYSTEM) $(LIBBOOST_IOSTREAMS)
test_barcode_SOURCES = test_barcode.cpp ../src/Barcode.h ../src/Barcode.cpp

test_demuxer_CPPFLAGS = -I../src/ ${BOOST_CPPFLAGS} $(BAMTOOLS_CFLAGS)
test_demuxer_LDADD = ${BOOST_LDFLAGS} ${LIBBOOST_UNIT_TEST_FRAMEWORK} $(LIBZ) $(LIBBOOST_FILESYSTEM) $(LIBBOOST_SYSTEM) $(LIBBOOST_IOSTREAMS) $(LIBDL) $(LIBBAMTOOLS) $(BAMTOOLS_LIBS)
test_demuxer_SOURCES = test_demuxer.cpp \
                       ../src/Barcode.h \
                       ../src/Barcode.cpp \
                       ../src/Parser.h \
                       ../src/Parser.cpp \
                       ../src/PairedReader.h \
                       ../src/PairedReader.cpp \
                       ../src/BoostZipReader.h \
                       ../src/BoostZipReader.cpp \
                       ../src/FastqReader.h \
                       ../src/FastqReader.cpp \
                       ../src/helper.h \
                       ../src/ZipFastqWriter.h \
                       ../src/ZipFastqWriter.cpp \
                       ../src/Writer.h \
                       ../src/Writer.cpp \
                       ../src/FileHandler.h \
                       ../src/FileHandler.cpp \
                       ../src/FileHandlerSE.h \
                       ../src/FileHandlerSE.cpp \
                       ../src/Correction_Counter.h \
                       ../src/Correction_Counter.cpp
if HAVE_LIBBAMTOOLS
    test_demuxer_SOURCES += ../src/BamReader.h \
                            ../src/BamReader.cpp
endif HAVE_LIBBAMTOOLS



test_parser_CPPFLAGS = -I../src/ ${BOOST_CPPFLAGS} $(BAMTOOLS_CFLAGS)
test_parser_LDADD = ${BOOST_LDFLAGS} ${LIBBOOST_UNIT_TEST_FRAMEWORK} $(LIBZ) $(LIBBOOST_FILESYSTEM) $(LIBBOOST_SYSTEM) $(LIBBOOST_IOSTREAMS) $(LIBDL) $(LIBBAMTOOLS) $(BAMTOOLS_LIBS)
test_parser_SOURCES = test_parser.cpp \
                      ../src/helper.h \
                      ../src/Barcode.h \
                      ../src/Barcode.cpp \
                      ../src/Parser.h \
                      ../src/Parser.cpp \
                      ../src/PairedReader.h \
                      ../src/PairedReader.cpp \
                      ../src/BoostZipReader.h \
                      ../src/BoostZipReader.cpp \
                      ../src/FastqReader.h \
                      ../src/FastqReader.cpp \
                      ../src/Correction_Counter.h \
                      ../src/Correction_Counter.cpp
if HAVE_LIBBAMTOOLS
    test_parser_SOURCES += ../src/BamReader.h \
                           ../src/BamReader.cpp
endif HAVE_LIBBAMTOOLS

