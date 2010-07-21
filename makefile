# Makefile for align
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
ZEEP_DIR			= $(HOME)/projects/libzeep/
MRS_LIB_DIR			= $(HOME)/projects/mrs/lib

DESTDIR				?= /usr/local/
LIBDIR				= $(DESTDIR)lib
INCDIR				= $(DESTDIR)include
MANDIR				= $(DESTDIR)man/man3

BOOST_LIBS			= system thread regex filesystem program_options
BOOST_LIBS			:= $(BOOST_LIBS:%=boost_%$(BOOST_LIB_SUFFIX))
LIBS				= zeep mrs $(BOOST_LIBS) z bz2 uuid 
LDOPTS				= $(BOOST_LIB_DIR:%=-L%) -L$(MRS_LIB_DIR) -L$(ZEEP_DIR) $(LIBS:%=-l%) -gdwarf-2 -pthread

CC					?= c++
CFLAGS				= $(BOOST_INC_DIR:%=-I%) -I$(ZEEP_DIR) -I$(MRS_LIB_DIR)/Sources \
					  -iquote ./ -gdwarf-2 -fPIC -O3 -pthread -Wno-multichar

VPATH += src

OBJECTS = \
	obj/align.o

align: $(OBJECTS)
	c++ -o $@ $(OBJECTS) $(LDOPTS)

obj/%.o: %.cpp
	c++ -MD -c -o $@ $< $(CFLAGS)

include $(OBJECTS:%.o=%.d)

$(OBJECTS:.o=.d):

clean:
	rm -rf obj/* align
