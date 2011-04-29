# Makefile for mas
#
#  Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
# Distributed under the Boost Software License, Version 1.0.
#    (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
#
# You may have to edit the first three defines on top of this
# makefile to match your current installation.

#BOOST_LIB_SUFFIX	= -mt
BOOST_LIB_DIR		= $(HOME)/projects/boost/lib
BOOST_INC_DIR		= $(HOME)/projects/boost/include
ZEEP_DIR			= $(HOME)/projects/libzeep/

DEST_DIR			?= /usr/local/
LIB_DIR				= $(DEST_DIR)lib $(HOME)/projects/mrs/lib \
					  $(ZEEP_DIR) $(BOOST_LIB_DIR)
INC_DIR				= $(BOOST_INC_DIR) $(HOME)/projects/mrs/lib/Sources \
					  $(ZEEP_DIR) src/
MAN_DIR				= $(DEST_DIR)man/man3

BOOST_LIBS			= thread regex filesystem program_options date_time iostreams math_c99 system
BOOST_LIBS			:= $(BOOST_LIBS:%=boost_%$(BOOST_LIB_SUFFIX))
LIBS				= zeep $(BOOST_LIBS) z bz2
LDOPTS				= $(LIB_DIR:%=-L%)
LDOPTS				+= $(LIBS:%=-l%) -gdwarf-2 -pthread

CC					?= c++
CFLAGS				= $(INC_DIR:%=-I%) -I$(ZEEP_DIR) -I$(MRS_LIB_DIR)/Sources -DBOOST_FILESYSTEM_VERSION=2 \
					  -iquote ./ -gdwarf-2 -Wall -Wno-multichar -pthread -std=c++0x
OPT					= -O3 -DNDEBUG # -march=native

include make.config

CFLAGS				+= $(OPT) -g -DLINUX -DUSE_COMPRESSION

VPATH += src

OBJECTS = \
	obj/align-3d.o \
	obj/dssp.o \
	obj/guide.o \
	obj/ioseq.o \
	obj/mas.o \
	obj/matrix.o \
	obj/primitives-3d.o \
	obj/structure.o \
	obj/utils.o

mas: $(OBJECTS)
	@ echo linking $@
	@ c++ -o $@ $(OBJECTS) $(LDOPTS)
	@ echo OK

dssp-2: obj/mkdssp.o obj/dssp.o obj/matrix.o obj/primitives-3d.o obj/structure.o obj/utils.o
	@ echo linking $@
	@ c++ -static -o $@ $? $(LDOPTS)
	@ echo OK

mkhssp: obj/mkhssp.o obj/dssp.o obj/matrix.o obj/primitives-3d.o obj/structure.o obj/utils.o
	@ echo linking $@
	@ c++ -o $@ $^ $(LDOPTS)
	@ echo OK

hsspsoap: obj/blast.o obj/dssp.o obj/hsspsoap.o obj/matrix.o obj/maxhom-hssp.o obj/primitives-3d.o obj/structure.o obj/utils.o
	@ echo linking $@
	@ c++ -o $@ $^ $(LDOPTS)
	@ echo OK

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
	sudo install -m 755 mas $(DEST_DIR)bin/mas

make.config:
	@echo "creating empty make.config file"
	@echo "# Set local options for make here" > make.config

