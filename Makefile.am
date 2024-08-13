AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS = -I m4

if CODE_COVERAGE_ENABLED
# top_builddir= .
CODE_COVERAGE_DIRECTORY=.
CODE_COVERAGE_BRANCH_COVERAGE=0
CODE_COVERAGE_LCOV_SHOPTS_DEFAULT=--no-external
# --no-recursion
include $(top_builddir)/aminclude_static.am
endif

# the subdirectories of the project to go into
SUBDIRS = src doc misc tests
pkgconfigdir = $(libdir)/pkgconfig

EXTRA_DIST =  packaging/debian \
              doc/man/man1/idemuxCPP.1 \
              README.md \
              LICENCE \
              misc/barcodes/* \
              tests/resources/*

if CODE_COVERAGE_ENABLED
clean-local: code-coverage-clean
distclean-local: code-coverage-dist-clean
endif
