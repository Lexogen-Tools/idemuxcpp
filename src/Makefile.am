AUTOMAKE_OPTIONS = subdir-objects

# the program to build
bin_PROGRAMS = 	idemuxCPP

SUFFIXES = _cmdl.c _cmdl.h .ggo

GENGETOPT_CMDL =  idemuxCPP_cmdl.c idemuxCPP_cmdl.h
GENGETOPT_FILES =  idemuxCPP.ggo

EXTRA_DIST =  ${GENGETOPT_FILES} ${GENGETOPT_CMDL}

CPPADD = @LIBZ@
LDADD = $(LIBZ) $(LIBBOOST_SYSTEM) $(LIBBOOST_FILESYSTEM)

# remove all generated content
maintainer-clean-local:
	rm -rf $(GENGETOPT_CMDL)

# how we make the commandline interface parser
.ggo_cmdl.c:
	$(AM_V_GEN)gengetopt -i $<

_cmdl.c_cmdl.h :

# idemuxCPP sources
idemuxCPP_SOURCES =  ${GENGETOPT_CMDL} \
                   idemuxCPP_cmdl.h \
                   idemuxCPP.cpp \
                   ZipFastqWriter.h \
                   ZipFastqWriter.cpp \
                   ZipFastqReader.h \
                   ZipFastqReader.cpp \
                   Writer.h \
                   Writer.cpp \
                   Parser.h \
                   Parser.cpp \
                   PairedReader.h \
                   PairedReader.cpp \
                   FileHandler.h \
                   FileHandler.cpp \
                   Barcode.h \
                   Barcode.cpp \
                   helper.h \
                   Demux.h




