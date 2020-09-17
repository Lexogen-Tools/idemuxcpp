AUTOMAKE_OPTIONS = subdir-objects
TESTS = test_barcode test_demuxer test_parser
check_PROGRAMS = test_barcode test_demuxer test_parser
test_barcode_CPPFLAGS = -I../src/ ${BOOST_CPPFLAGS}
test_barcode_LDADD = ${BOOST_LDFLAGS} ${LIBBOOST_UNIT_TEST_FRAMEWORK} $(LIBBOOST_FILESYSTEM_LIB) $(LIBBOOST_SYSTEM)
test_barcode_SOURCES = test_barcode.cpp ../src/Barcode.h ../src/Barcode.cpp

test_demuxer_CPPFLAGS = -I../src/ ${BOOST_CPPFLAGS}
test_demuxer_LDADD = ${BOOST_LDFLAGS} ${LIBBOOST_UNIT_TEST_FRAMEWORK} $(LIBZ) $(LIBBOOST_FILESYSTEM) $(LIBBOOST_SYSTEM) $(LIBDL)
test_demuxer_SOURCES = test_demuxer.cpp \
                       ../src/Barcode.h \
                       ../src/Barcode.cpp \
                       ../src/Parser.h \
                       ../src/Parser.cpp \
                       ../src/PairedReader.h \
                       ../src/PairedReader.cpp \
                       ../src/ZipFastqReader.h \
                       ../src/ZipFastqReader.cpp \
                       ../src/helper.h \
                       ../src/ZipFastqWriter.h \
                       ../src/ZipFastqWriter.cpp \
                       ../src/Writer.h \
                       ../src/Writer.cpp \
                       ../src/FileHandler.h \
                       ../src/FileHandler.cpp

test_parser_CPPFLAGS = -I../src/ ${BOOST_CPPFLAGS}
test_parser_LDADD = ${BOOST_LDFLAGS} ${LIBBOOST_UNIT_TEST_FRAMEWORK} $(LIBZ) $(LIBBOOST_FILESYSTEM) $(LIBBOOST_SYSTEM) $(LIBDL)
test_parser_SOURCES = test_parser.cpp \
                      ../src/helper.h \
                      ../src/Barcode.h \
                      ../src/Barcode.cpp \
                      ../src/Parser.h \
                      ../src/Parser.cpp \
                      ../src/PairedReader.h \
                      ../src/PairedReader.cpp \
                      ../src/ZipFastqReader.h \
                      ../src/ZipFastqReader.cpp


