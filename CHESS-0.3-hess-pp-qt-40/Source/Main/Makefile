# This file is part of CHESS.
#      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
#      2013
#
# The file is modified from Makefile in the ESS++ program
#      Copyright (c) Marc Chadeau (m.chadeau@imperial.ac.uk)
#                    Leonardo Bottolo (l.bottolo@imperial.ac.uk)
#                    David Hastie (d.hastie@imperial.ac.uk)
#      2010
#
# CHESS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CHESS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CHESS.  If not, see <http://www.gnu.org/licenses/>.


CUDA			= 0
GOTO                    = 0

LIBDIR			= ../Kernel/Routines/
LIBDIR2			= ../Kernel/Classes/
LIBDIR3			= ../General_Classes/
LIBDIR4			= ../Models/ESS/
LIBDIR5			= ../Models/HESS/

ifeq ($(GOTO),1)
        #GOTODIR		= /path-to-GOTO-library
endif


CC	 		= gcc
CXX 		= g++

IFLAGSTMP		= -I$(LIBDIR)

ifeq ($(GOTO),1)
		LFLAGSTMP		= -L$(LIBDIR) -L$(GOTODIR)
else
		LFLAGSTMP		= -L$(LIBDIR)
endif		

LIBSTMP			= -lgsl

ifeq ($(CUDA),0)
	IFLAGS		= $(IFLAGSTMP)
	LFLAGS 		= $(LFLAGSTMP)
        ifeq ($(GOTO),1)
            LIBS 		= $(LIBSTMP) -lgoto2
        else
            LIBS 		= $(LIBSTMP) -lgslcblas
        endif
else
	IFLAGS		= $(IFLAGSTMP) -I/usr/local/cuda/include -I/usr/local/cula/include
	LFLAGS		= $(LFLAGSTMP) -L/usr/local/cuda/lib64 -L/usr/local/cula/lib64
	LIBS		= $(LIBSTMP) -lgoto2 -lpthread -llapack -lcuda -lcublas -lcudart -lcula
endif

CXXFLAGS		= 
CPPFLAGS		= -W -Wall --pedantic -O3  $(IFLAGS) -D_CUDA_=$(CUDA)
#CPPFLAGS		= -W -Wall -g -pg -O0 --pedantic $(IFLAGS) -D_CUDA_=$(CUDA)
#CPPFLAGS		= -W -Wall -O4 -g $(IFLAGS) -D_CUDA_=$(CUDA)
# Use the compile option -pg for profiling.


all:		CHESS .depend 

clean:		
		rm -f *.o
		rm -f .depend
		rm -f CHESS
		rm -f *~
		rm -f ../Kernel/Routines/*.o
		rm -f ../Kernel/Routines/*~
		rm -f ../Kernel/Classes/*.o
		rm -f ../Kernel/Classes/*~
		rm -f ../General_Classes/*.o
		rm -f ../General_Classes/*~
		rm -f ../Models/ESS/*.o
		rm -f ../Models/ESS/*~
		rm -f ../Models/HESS/*.o
		rm -f ../Models/HESS/*~
		rm -f ../.depend
dep .depend:
	$(CXX) $(CPPFLAGS) -MM CHESS.cc >>.depend

infiles1	= CHESS.o  $(LIBDIR)dyn_name.o $(LIBDIR)matrix_handling.o  $(LIBDIR)post_processing.o $(LIBDIR)rand.o $(LIBDIR2)String_Matrices.o $(LIBDIR2)Int_Matrices.o $(LIBDIR2)Double_Matrices.o  $(LIBDIR2)Double_Matrices_cont.o  $(LIBDIR)regression.o $(LIBDIR)cond_post.o $(LIBDIR)xml_file_read.o $(LIBDIR2)Prior_param.o $(LIBDIR)moves.o $(LIBDIR2)Temperatures.o $(LIBDIR2)AdMH.o $(LIBDIR2)Move_monitor.o $(LIBDIR2)DR.o $(LIBDIR2)CM.o $(LIBDIR2)Int_Matrices_var_dim.o $(LIBDIR3)Kernel_Single_Gamma.o $(LIBDIR3)Model_Generic.o $(LIBDIR4)Model_ESS.o $(LIBDIR4)Start_Up_ESS.o $(LIBDIR3)Model_Information.o $(LIBDIR5)Model_HESS.o $(LIBDIR3)Command_Line.o

CHESS:	$(infiles1)
		$(CXX) $(CPPFLAGS) $(LFLAGS) -o CHESS $(infiles1) $(LIBS)

