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

VERSION				= 2.0.3

#BOOST_LIB_SUFFIX	= -mt
BOOST_LIB_DIR		?= $(HOME)/projects/boost/lib
BOOST_INC_DIR		?= $(HOME)/projects/boost/include
ZEEP_DIR			?= $(HOME)/projects/libzeep/

DEST_DIR			?= /usr/local/
LIB_DIR				= $(BOOST_LIB_DIR) $(ZEEP_DIR) $(DEST_DIR)lib $(HOME)/projects/mrs/lib
INC_DIR				= $(BOOST_INC_DIR) $(HOME)/projects/mrs/lib/Sources \
					  $(ZEEP_DIR) src/
MAN_DIR				= $(DEST_DIR)man/man3

BOOST_LIBS			= thread regex filesystem program_options date_time iostreams math_c99 system
BOOST_LIBS			:= $(BOOST_LIBS:%=boost_%$(BOOST_LIB_SUFFIX))
LIBS				= zeep $(BOOST_LIBS) z bz2
ifeq ($(DEBUG),1)
LIBS				+= mrsd
else
LIBS				+= mrs
endif

LDOPTS				= $(LIB_DIR:%=-L%)
LDOPTS				+= $(LIBS:%=-l%) -gdwarf-2 -pthread

CC					= c++
CFLAGS				= $(INC_DIR:%=-I%) -I$(ZEEP_DIR) -I$(MRS_LIB_DIR)/Sources -DBOOST_FILESYSTEM_VERSION=2 \
					  -iquote ./ -gdwarf-2 -Wall -Wno-multichar -pthread \
					  -DBOOST_FILESYSTEM_VERSION=2 -std=c++0x -DVERSION='"$(VERSION)"'

ifneq ($(DEBUG),1)
OPT					= -O3 -DNDEBUG # -march=native
endif

include make.config

CFLAGS				+= $(OPT) -g -DLINUX -DUSE_COMPRESSION

OBJ					= obj

ifeq ($(DEBUG),1)
CFLAGS				+= -g3
OBJ					:= $(OBJ).dbg
endif

VPATH += src

SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:src/%.cpp=$(OBJ)/%.o)

all: mkdssp mkhssp sto2fa # hsspsoap

mas: $(OBJECTS)
	@ echo linking $@
	@ $(CC) -o $@ $(OBJECTS) $(LDOPTS)
	@ echo OK

mkdssp: $(OBJ)/mkdssp.o $(OBJ)/dssp.o $(OBJ)/primitives-3d.o $(OBJ)/structure.o $(OBJ)/utils.o
	@ echo linking $@
	@ $(CC) -static -o $@ $^ $(LDOPTS)
	@ echo OK

mkhssp: $(OBJ)/mkhssp.o $(OBJ)/dssp.o $(OBJ)/hmmer-hssp.o $(OBJ)/matrix.o $(OBJ)/primitives-3d.o $(OBJ)/structure.o $(OBJ)/utils.o
	@ echo linking $@
	@ $(CC) -o $@ $^ $(LDOPTS)
	@ echo OK

sto2fa: $(OBJ)/sto2fa.o $(OBJ)/dssp.o $(OBJ)/hmmer-hssp.o $(OBJ)/matrix.o $(OBJ)/primitives-3d.o $(OBJ)/structure.o $(OBJ)/utils.o
	@ echo linking $@
	@ $(CC) -o $@ $^ $(LDOPTS)
	@ echo OK

hsspsoap: $(OBJ)/dssp.o $(OBJ)/hsspsoap.o $(OBJ)/matrix.o $(OBJ)/maxhom-hssp.o $(OBJ)/primitives-3d.o $(OBJ)/structure.o $(OBJ)/utils.o
	@ echo linking $@
	@ $(CC) -o $@ $^ $(LDOPTS)
	@ echo OK

$(OBJ)/%.o: %.cpp
	@ if [ ! -d $(OBJ) ]; then mkdir $(OBJ); fi
	@ echo compiling $@
	@ $(CC) -MD -c -o $@ $< $(CFLAGS)

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

test:
	echo $(OBJECTS)
