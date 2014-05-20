bin_PROGRAMS	=	mkdssp mkhssp hsspconv hsspsoap

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
									src/iocif.cpp \
									src/mas.cpp \
									src/mkdssp.cpp \
									src/primitives-3d.cpp \
									src/structure.cpp \
									src/utils.cpp

mkdssp_LDADD	=	$(shared_LDADD)

mkhssp_SOURCES =	src/blast.cpp \
									src/dssp.cpp \
									src/fetchdbrefs.cpp \
									src/hssp-nt.cpp \
									src/iocif.cpp \
									src/mas.cpp \
									src/matrix.cpp \
									src/mkhssp.cpp \
									src/primitives-3d.cpp \
									src/progress.cpp \
									src/structure.cpp \
									src/utils.cpp

mkhssp_LDADD	=	$(shared_LDADD) \
								-lboost_timer

hsspconv_SOURCES = src/hssp-convert-3to1.cpp \
									 src/utils.cpp

hsspconv_LDADD	=	$(shared_LDADD)

hsspsoap_SOURCES =	src/blast.cpp \
										src/dssp.cpp \
										src/fetchdbrefs.cpp \
										src/hssp-nt.cpp \
										src/hsspsoap.cpp \
										src/iocif.cpp \
										src/mas.cpp \
										src/matrix.cpp \
										src/primitives-3d.cpp \
										src/progress.cpp \
										src/structure.cpp \
										src/utils.cpp

hsspsoap_LDADD	=	$(shared_LDADD) \
									-lboost_math_c99 \
									-lbz2 -lz -lrt \
									-lboost_timer \
									rsrc/index.o \
									rsrc/error.o

AM_CPPFLAGS	=	-std=c++11 \
							-pedantic \
							-Wall \
							-Werror \
							-Wno-reorder \
							$(BOOST_CPPFLAGS) \
							-Isrc/

BUILT_SOURCES	=	mtrx/matrices.h

rsrc/index.o: rsrc/index.html
	ld -r -b binary -o rsrc/index.o rsrc/index.html

rsrc/error.o: rsrc/error.html
	ld -r -b binary -o rsrc/error.o rsrc/error.html

mtrx/matrices.h: mtrx/mkmat_h.pl
	perl mtrx/mkmat_h.pl mtrx/

if HAVE_CPPCHECK
cppcheck:
		$(CPPCHECK) --quiet --enable=style src test
check:	cppcheck
endif
