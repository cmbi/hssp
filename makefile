# Makefile for the DSSP/HSSP software suite
#
#  Copyright Maarten L. Hekkelman, Radboud University 2008-2011.
# Distributed under the Boost Software License, Version 1.0.
#    (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)
#
# You may have to edit the first three defines on top of this
# makefile to match your current installation.

firstTarget: all

include make.config

VERSION				= 2.1.0

DEST_DIR			?= /usr/local/
LIB_DIR				= $(BOOST_LIB_DIR) $(ZEEP_DIR) $(DEST_DIR)lib $(HOME)/projects/mrs/lib
INC_DIR				= $(BOOST_INC_DIR) $(HOME)/projects/mrs/lib/Sources \
					  $(ZEEP_DIR) src/
MAN_DIR				= $(DEST_DIR)man/man3

BOOST_LIBS			= thread regex filesystem program_options date_time iostreams math_c99 system 
BOOST_LIBS			:= $(BOOST_LIBS:%=boost_%$(BOOST_LIB_SUFFIX))
LIBS				= zeep $(BOOST_LIBS) z bz2 rt
#ifeq ($(DEBUG),1)
#LIBS				+= mrsd
#else
#LIBS				+= mrs
#endif

LDOPTS				= $(LIB_DIR:%=-L%)
LDOPTS				+= $(LIBS:%=-l%) -gdwarf-2 -pthread

CC					?= icpc
CFLAGS				= $(INC_DIR:%=-I%) -I$(ZEEP_DIR) -I$(MRS_LIB_DIR)/Sources \
					  -iquote ./ -gdwarf-2 -Wall -Wno-multichar -pthread \
					  -std=c++0x -DVERSION='"$(VERSION)"'

ifneq ($(DEBUG),1)
OPT					= -O3 -DNDEBUG # -march=native
endif


CFLAGS				+= $(OPT) -g -DLINUX -DUSE_COMPRESSION

OBJ					= obj

ifeq ($(DEBUG),1)
CFLAGS				+= -g3
OBJ					:= $(OBJ).dbg
endif

VPATH += src $(OBJ)

SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:src/%.cpp=$(OBJ)/%.o)

all: mkdssp mkhssp sto2fa aln2hssp  | $(OBJ)

mas: $(OBJECTS)
	@ echo linking $@
	@ $(CC) -o $@ $(OBJECTS) $(LDOPTS)
	@ echo OK

mkdssp: mkdssp.o dssp.o primitives-3d.o structure.o utils.o mas.o
	@ echo linking $@
	@ $(CC) -static -o $@ $^ $(LDOPTS)
	@ echo OK

mkhssp: mkhssp.o hssp-nt.o dssp.o matrix.o primitives-3d.o structure.o utils.o blast.o mas.o matrix.o fetchdbrefs.o
	@ echo linking $@
	@ $(CC) $(CFLAGS) -o $@ $^ $(LDOPTS) 
	@ echo OK

hsspsoap: hsspsoap.o hssp-nt.o dssp.o matrix.o primitives-3d.o structure.o utils.o blast.o mas.o matrix.o fetchdbrefs.o
	@ echo linking $@
	@ $(CC) -static -o $@ $^ $(LDOPTS)
	@ echo OK

hsspconv: hssp-convert-3to1.o utils.o
	@ echo linking $@
	@ $(CC) -static -o $@ $^ $(LDOPTS)
	@ echo OK

$(OBJ)/%.o: %.cpp | $(OBJ)
	@ if [ ! -d $(OBJ) ]; then mkdir $(OBJ); fi
	@ echo compiling $@
	@ $(CC) -MD -c -o $@ $< $(CFLAGS)

$(OBJ):
	mkdir -p $@

$(OBJ)/matrix.o: mtrx/matrices.h

mtrx/matrices.h: mtrx/mkmat_h.pl
	perl mtrx/mkmat_h.pl mtrx/

include $(OBJECTS:%.o=%.d)

$(OBJECTS:.o=.d):

clean:
	rm -rf $(OBJ)/* mas

install: mas
	sudo install -m 755 mas $(DEST_DIR)bin/mas

make.config:
	@echo "creating empty make.config file"
	@echo "# Set local options for make here" > make.config
	@echo "#BOOST_LIB_SUFFIX	= -mt" >> make.config
	@echo "#BOOST_LIB_DIR		?= $(HOME)/projects/boost/lib" >> make.config
	@echo "#BOOST_INC_DIR		?= $(HOME)/projects/boost/include" >> make.config
	@echo "#ZEEP_DIR			?= $(HOME)/projects/libzeep/" >> make.config


test:
	echo $(OBJECTS)
