# Makefile for mas
#
#  Copyright Maarten L. Hekkelman, Radboud University 2008-2010.
# Distributed under the Boost Software License, Version 1.0.
#    (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
#
# You may have to edit the first three defines on top of this
# makefile to match your current installation.

#BOOST_LIB_SUFFIX	= -mt				# Works for Ubuntu
BOOST_LIB_DIR		= $(HOME)/projects/boost/lib
BOOST_INC_DIR		= $(HOME)/projects/boost/include
#ZEEP_DIR			= $(HOME)/projects/libzeep/
#MRS_LIB_DIR			= $(HOME)/projects/mrs/lib

DESTDIR				?= /usr/local/
LIBDIR				= $(DESTDIR)lib
INCDIR				= $(DESTDIR)include
MANDIR				= $(DESTDIR)man/man3

BOOST_LIBS			= system thread regex filesystem program_options
BOOST_LIBS			:= $(BOOST_LIBS:%=boost_%$(BOOST_LIB_SUFFIX))
#LIBS				= zeep mrs $(BOOST_LIBS) z bz2 uuid 
LIBS				= $(BOOST_LIBS)
LDOPTS				= $(BOOST_LIB_DIR:%=-L%) # -L$(MRS_LIB_DIR) -L$(ZEEP_DIR)
LDOPTS				+= $(LIBS:%=-l%) -gdwarf-2 -pthread

CC					?= c++
CFLAGS				= $(BOOST_INC_DIR:%=-I%) -I$(ZEEP_DIR) -I$(MRS_LIB_DIR)/Sources \
					  -iquote ./ -gdwarf-2 -fPIC -pthread -Wno-multichar -std=c++0x
#CFLAGS				+= -O3  # -DNDEBUG

VPATH += src

OBJECTS = \
	obj/ioseq.o \
	obj/mas.o \
	obj/matrix.o \
	obj/utils.o

mas: $(OBJECTS)
	@ echo linking $@
	@ c++ -o $@ $(OBJECTS) $(LDOPTS)

obj/%.o: %.cpp
	@ echo compiling $@
	@ c++ -MD -c -o $@ $< $(CFLAGS)

obj/matrix.o: mtrx/matrices.h

mtrx/matrices.h: mtrx/mkmat_h.pl
	perl mtrx/mkmat_h.pl mtrx/

include $(OBJECTS:%.o=%.d)

$(OBJECTS:.o=.d):

clean:
	rm -rf obj/* mas

install: mas
	sudo install -m 755 mas $(DESTDIR)bin/mas

