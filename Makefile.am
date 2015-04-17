bin_PROGRAMS	=	mkdssp mkhssp hsspconv hsspsoap test_xssp

shared_LDADD =	-lzeep \
								$(BOOST_DATE_TIME_LIB) \
								$(BOOST_FILESYSTEM_LIB) \
								$(BOOST_IOSTREAMS_LIB) \
								$(BOOST_PROGRAM_OPTIONS_LIB) \
								$(BOOST_REGEX_LIB) \
								$(BOOST_SYSTEM_LIB) \
								$(BOOST_THREAD_LIB) \
								-lpthread

mkdssp_SOURCES	=	src/dssp.cpp \
									src/dssp.h \
									src/iocif.cpp \
									src/iocif.h \
									src/mas.cpp \
									src/mas.h \
									src/mkdssp.cpp \
									src/primitives-3d.cpp \
									src/primitives-3d.h \
									src/structure.cpp \
									src/structure.h \
									src/buffer.h \
									src/utils.cpp \
									src/utils.h \
									src/version.h

mkdssp_LDADD	=	$(shared_LDADD)

mkhssp_SOURCES =	src/blast.cpp \
									src/blast.h \
									src/dssp.cpp \
									src/dssp.h \
									src/fasta.cpp \
									src/fasta.h \
									src/fetchdbrefs.cpp \
									src/fetchdbrefs.h \
									src/hssp-nt.cpp \
									src/hssp-nt.h \
									src/iocif.cpp \
									src/iocif.h \
									src/mas.cpp \
									src/mas.h \
									src/matrix.cpp \
									src/matrix.h \
									mtrx/matrices.h \
									mtrx/mkmat_h.pl \
									src/mkhssp.cpp \
									src/primitives-3d.cpp \
									src/primitives-3d.h \
									src/progress.cpp \
									src/progress.h \
									src/structure.cpp \
									src/structure.h \
									src/buffer.h \
									src/align-2d.h \
									src/utils.cpp \
									src/utils.h

mkhssp_LDADD	=	$(shared_LDADD) \
								-lboost_timer

hsspconv_SOURCES = src/hssp-convert-3to1.cpp \
									 src/utils.cpp \
									 src/utils.h


hsspconv_LDADD	=	$(shared_LDADD)

hsspsoap_SOURCES =	src/blast.cpp \
										src/blast.h \
										src/dssp.cpp \
										src/dssp.h \
										src/fetchdbrefs.cpp \
										src/fetchdbrefs.h \
										src/hssp-nt.cpp \
										src/hssp-nt.h \
										src/hsspsoap.cpp \
										src/iocif.cpp \
										src/iocif.h \
										src/mas.cpp \
										src/mas.h \
										src/matrix.cpp \
										src/matrix.h \
										src/primitives-3d.cpp \
										src/primitives-3d.h \
										src/progress.cpp \
										src/progress.h \
										src/structure.cpp \
										src/structure.h \
										src/buffer.h \
										src/align-2d.h \
										src/utils.cpp \
										src/utils.h \
										rsrc/index.html \
										rsrc/error.html

hsspsoap_LDADD	=	$(shared_LDADD) \
									-lboost_math_c99 \
									-lbz2 -lz -lrt \
									-lboost_timer \
									rsrc/index.o \
									rsrc/error.o

test_xssp_SOURCES	=	src/fasta.cpp \
										src/iocif.cpp \
										src/mas.cpp \
										src/primitives-3d.cpp \
										src/primitives-3d.h \
										src/structure.cpp \
										src/structure.h \
										src/buffer.h \
										src/align-2d.h \
										src/utils.cpp \
										tests/test_fasta.cpp

test_xssp_LDADD	=	$(shared_LDADD) \
									$(BOOST_UNIT_TEST_FRAMEWORK_LIB)

AM_CPPFLAGS	=	-std=c++0x \
							-pedantic \
							-Wall \
							-Werror \
							-Wno-reorder \
							$(BOOST_CPPFLAGS) \
							-Isrc/

BUILT_SOURCES	=	mtrx/matrices.h
CLEANFILES = mtrx/matrices.h

rsrc/index.o: rsrc/index.html
	ld -r -b binary -o rsrc/index.o rsrc/index.html

rsrc/error.o: rsrc/error.html
	ld -r -b binary -o rsrc/error.o rsrc/error.html

mtrx/matrices.h: mtrx/mkmat_h.pl
	$(PERL) mtrx/mkmat_h.pl mtrx/

if HAVE_CPPCHECK
cppcheck:
		$(CPPCHECK) --quiet --enable=style src test
check:	cppcheck
endif
