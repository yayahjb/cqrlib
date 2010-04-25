#
#  Makefile
#  CQRlib
#
#  Created by Herbert J. Bernstein on 02/22/09.
#  Copyright 2009 Herbert J. Bernstein. All rights reserved.
#
#

#  Work supported in part by NIH NIGMS under grant 1R15GM078077-01 and DOE 
#  under grant ER63601-1021466-0009501.  Any opinions, findings, and 
#  conclusions or recommendations expressed in this material are those of the   
#  author(s) and do not necessarily reflect the views of the funding agencies.


#**********************************************************************
#*                                                                    *
#* YOU MAY REDISTRIBUTE THE CQRlib API UNDER THE TERMS OF THE LGPL    *
#*                                                                    *
#**********************************************************************/

#************************* LGPL NOTICES *******************************
#*                                                                    *
#* This library is free software; you can redistribute it and/or      *
#* modify it under the terms of the GNU Lesser General Public         *
#* License as published by the Free Software Foundation; either       *
#* version 2.1 of the License, or (at your option) any later version. *
#*                                                                    *
#* This library is distributed in the hope that it will be useful,    *
#* but WITHOUT ANY WARRANTY; without even the implied warranty of     *
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  *
#* Lesser General Public License for more details.                    *
#*                                                                    *
#* You should have received a copy of the GNU Lesser General Public   *
#* License along with this library; if not, write to the Free         *
#* Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,    *
#* MA  02110-1301  USA                                                *
#*                                                                    *
#**********************************************************************/

# Version string
VERSION = 2:0:0
RELEASE = 1.0.3


#
# Compiler and compilation flags
#
CC	= gcc
CXX = g++
CFLAGS  = -g -O2  -Wall -ansi -pedantic
CPPFLAGS = $(CFLAGS) -DCQR_NOCCODE=1

#
# libtool path if system default is not suitable
#
#LIBTOOL = $(HOME)/bin/libtool
ifndef LIBTOOL
  LIBTOOL = libtool
endif

#
# If local headers must be quoted uncomment the next line
#
#USE_LOCAL_HEADERS = 1


#
# Directories
#
ROOT     = .
LIB      = $(ROOT)/lib
BIN      = $(ROOT)/bin
SRC      = $(ROOT)
INC      = $(ROOT)
EXAMPLES = $(ROOT)
TESTDATA = $(ROOT)
#INSTALLDIR = /usr/local
INSTALLDIR  = $(HOME)

#
# Include directories
#
ifdef USE_LOCAL_HEADERS
INCLUDES = -DUSE_LOCAL_HEADERS
else
INCLUDES = -I$(INC)
endif

COMPILE_COMMAND         =  $(LIBTOOL) --mode=compile $(CC) $(CFLAGS) $(INCLUDES) $(WARNINGS) -c
LIBRARY_LINK_COMMAND    =  $(LIBTOOL) --mode=link  $(CC) -version-info $(VERSION) -rpath $(INSTALLDIR)/lib
BUILD_COMMAND_LOCAL     =  $(LIBTOOL) --mode=link $(CC) $(CFLAGS) $(INCLUDES)
BUILD_COMMAND_DYNAMIC   =  $(LIBTOOL) --mode=link $(CC) $(CFLAGS) -dynamic -I $(INSTALLDIR)/include -L$(INSTALLDIR)/lib
BUILD_COMMAND_STATIC    =  $(LIBTOOL) --mode=link $(CC) $(CFLAGS) -static -I $(INSTALLDIR)/include -L$(INSTALLDIR)/lib
CPPCOMPILE_COMMAND      =  $(LIBTOOL) --mode=compile $(CXX) $(CPPFLAGS) $(INCLUDES) $(WARNINGS) -c
CPPLIBRARY_LINK_COMMAND =  $(LIBTOOL) --mode=link $(CXX) -version-info $(VERSION) -rpath $(INSTALLDIR)/lib
CPPBUILD_COMMAND_LOCAL  =  $(LIBTOOL) --mode=link $(CXX) $(CPPFLAGS) $(INCLUDES)
CPPBUILD_COMMAND_DYNAMIC=  $(LIBTOOL) --mode=link $(CXX) $(CPPFLAGS) -dynamic -I $(INSTALLDIR)/include -L$(INSTALLDIR)/lib
CPPBUILD_COMMAND_STATIC =  $(LIBTOOL) --mode=link $(CXX) $(CPPFLAGS) -static -I $(INSTALLDIR)/include -L$(INSTALLDIR)/lib
INSTALL_COMMAND         =  $(LIBTOOL) --mode=install cp
INSTALL_FINISH_COMMAND  =  $(LIBTOOL) --mode=finish

OBJ_EXT                =  lo
LIB_EXT                =  la

######################################################################
#  You should not need to make modifications below this line         #
######################################################################

#
# Suffixes of files to be used or built
#
.SUFFIXES:	.c .$(OBJ_EXT) .$(LIB_EXT)

#
# Common dependencies
#
COMMONDEP = Makefile

#
# Source files
#
SOURCE   =  $(SRC)/cqrlib.c

#
# Header files
#
HEADERS   =  $(INC)/cqrlib.h

# Default: instructions
#
default:
	@echo ' '
	@echo '***************************************************************'
	@echo ' '
	@echo ' PLEASE READ README_CQRlib.txt and lgpl.txt'
	@echo ' '
	@echo ' Before making the CQRlib library and example programs, check'
	@echo ' that the chosen settings are correct'
	@echo ' '
	@echo ' The current compile command is:'
	@echo ' '
	@echo '   $(COMPILE_COMMAND)'
	@echo ' '
	@echo ' The current library link command is:'
	@echo ' '
	@echo '   $(LIBRARY_LINK_COMMAND)'
	@echo ' '
	@echo ' The current library local, dynamic and static build commands are:'
	@echo ' '
	@echo '   $(BUILD_COMMAND_LOCAL)'
	@echo '   $(BUILD_COMMAND_DYNAMIC)'
	@echo '   $(BUILD_COMMAND_STATIC)'
	@echo ' '
	@echo ' Before installing the CQRlib library and example programs, check'
	@echo ' that the install directory and install commands are correct:'
	@echo ' '
	@echo ' The current values are :'
	@echo ' '
	@echo '   $(INSTALLDIR) '	
	@echo '   $(INSTALL_COMMAND) '	
	@echo '   $(INSTALL_FINISH) '	
	@echo ' '
	@echo ' To compile the CQRlib library and example programs type:'
	@echo ' '
	@echo '   make clean'
	@echo '   make all'
	@echo ' '
	@echo ' To run a set of tests type:'
	@echo ' '
	@echo '   make tests'
	@echo ' '
	@echo ' To clean up the directories type:'
	@echo ' '
	@echo '   make clean'
	@echo ' '
	@echo ' To install the library and binaries type:'
	@echo ' '
	@echo '   make install'
	@echo ' '
	@echo '***************************************************************'
	@echo ' '

#
# Compile the library and examples
#
all:	$(LIB) $(BIN) $(SOURCE) $(HEADERS) \
		$(LIB)/libCQRlib.$(LIB_EXT) \
		$(BIN)/CQRlibTest $(BIN)/CPPQRTest
		
install:  all $(INSTALLDIR) $(INSTALLDIR)/lib $(INSTALLDIR)/include \
          $(INC) $(LIB)/libCQRlib.$(LIB_EXT)  $(INC)/cqrlib.h
		  $(INSTALL_COMMAND) $(LIB)/libCQRlib.$(LIB_EXT) $(INSTALLDIR)/lib/libCQRlib.$(LIB_EXT)
		  $(INSTALL_FINISH_COMMAND) $(INSTALLDIR)/lib/libCQRlib.$(LIB_EXT)
		  -cp $(INSTALLDIR)/include/cqrlib.h $(INSTALLDIR)/include/CQRlib_old.h
		  cp $(INC)/cqrlib.h $(INSTALLDIR)/include/cqrlib.h
		  chmod 644 $(INSTALLDIR)/include/cqrlib.h
		  echo "Testing final install dynamic"
		  $(BUILD_COMMAND_DYNAMIC)  $(EXAMPLES)/CQRlibTest.c \
		  -lCQRlib -lm -o $(BIN)/CQRlibTest_dynamic
		  $(BIN)/CQRlibTest_dynamic > $(TESTDATA)/CQRlibTest_dynamic.lst
		  diff -b -c $(TESTDATA)/CQRlibTest_orig.lst \
		    $(TESTDATA)/CQRlibTest_dynamic.lst
		  echo "Testing final install static"
		  $(BUILD_COMMAND_STATIC) $(EXAMPLES)/CQRlibTest.c \
		   -lCQRlib -lm -o $(BIN)/CQRlibTest_static
		  $(BIN)/CQRlibTest_static > $(TESTDATA)/CQRlibTest_static.lst
		  diff -b -c $(TESTDATA)/CQRlibTest_orig.lst \
		    $(TESTDATA)/CQRlibTest_static.lst
		  $(CPPBUILD_COMMAND_DYNAMIC) $(EXAMPLES)/CPPQRTest.cpp \
		    -lm -o $(BIN)/CPPQRTest_dynamic
		  $(BIN)/CPPQRTest_dynamic > $(TESTDATA)/CPPQRTest_dynamic.lst
		  diff -b -c $(TESTDATA)/CPPQRTest_orig.lst \
		    $(TESTDATA)/CPPQRTest_dynamic.lst
		  $(CPPBUILD_COMMAND_STATIC) $(EXAMPLES)/CPPQRTest.cpp \
		    -lm -o $(BIN)/CPPQRTest_static
	      $(BIN)/CPPQRTest_static > $(TESTDATA)/CPPQRTest_static.lst
		  diff -b -c $(TESTDATA)/CPPQRTest_orig.lst \
		    $(TESTDATA)/CPPQRTest_static.lst
			
		  
		  
#
# Directories
#
$(INSTALLDIR):
	mkdir -p $(INSTALLDIR)

$(INSTALLDIR)/lib:  $(INSTALLDIR)
	mkdir -p $(INSTALLDIR)/lib

$(INSTALLDIR)/bin:  $(INSTALLDIR)
	mkdir -p $(INSTALLDIR)/bin

$(INSTALLDIR)/include:  $(INSTALLDIR)
	mkdir -p $(INSTALLDIR)/include
	

$(LIB):
	mkdir $(LIB)

$(BIN):
	mkdir $(BIN)

#
# CQRlib library
#
$(LIB)/libCQRlib.$(LIB_EXT): $(SOURCE) $(HEADERS) $(COMMONDEP)
	$(COMPILE_COMMAND) -c $(SOURCE)
	$(LIBRARY_LINK_COMMAND) -o $(LIB)/libCQRlib.$(LIB_EXT) *.$(OBJ_EXT) 

#
# CQRlibTest example program
#
$(BIN)/CQRlibTest: $(LIB)/libCQRlib.$(LIB_EXT) $(EXAMPLES)/CQRlibTest.c 
	$(BUILD_COMMAND_LOCAL) $(EXAMPLES)/CQRlibTest.c $(LIB)/libCQRlib.$(LIB_EXT) -lm \
		      -o $@

#
# CPPQRTest example program
#
$(BIN)/CPPQRTest: $(EXAMPLES)/CPPQRTest.cpp
	$(CPPBUILD_COMMAND_LOCAL) $(EXAMPLES)/CPPQRTest.cpp -lm \
		      -o $@

#
# Tests
#
tests:		$(LIB) $(BIN) $(BIN)/CQRlibTest \
				all $(TESTDATA)/CQRlibTest_orig.lst $(TESTDATA)/CPPQRTest_orig.lst
			$(BIN)/CQRlibTest > $(TESTDATA)/CQRlibTest.lst
			diff -b -c $(TESTDATA)/CQRlibTest_orig.lst \
				$(TESTDATA)/CQRlibTest.lst
			$(BIN)/CPPQRTest > $(TESTDATA)/CPPQRTest.lst
			diff -b -c $(TESTDATA)/CPPQRTest_orig.lst \
				$(TESTDATA)/CPPQRTest.lst
			

#
# Remove all non-source files
#
empty:
		  @-rm -rf $(LIB)
		  @-rm -rf $(BIN)
		  @-rm -f $(TESTDATA)/CQRlibTest.lst
		  @-rm -f $(TESTDATA)/CQRlibTest_static.lst
		  @-rm -f $(TESTDATA)/CQRlibTest_dynamic.lst
		  @-rm -f $(TESTDATA)/CPPQRTest.lst
		  @-rm -f $(TESTDATA)/CPPQRTest_static.lst
		  @-rm -f $(TESTDATA)/CPPQRTest_dynamic.lst
		  
#
# Remove temporary files
#
clean:
	@-rm -f core 
	@-rm -f *.o
	@-rm -f *.$(OBJ_EXT)
	@-rm -f *.c.*

#
# Restore to distribution state
#
distclean:	clean empty





