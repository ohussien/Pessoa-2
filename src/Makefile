#----------------------------------------------------------------------------#
# Makefile for Pessoa source files								     #
# University of California, Los Angeles                                      #
# Version 2, Novemeber, 2015	                                     	     #
#----------------------------------------------------------------------------#



  CXX        = g++
  CFLAGS     = -pg -fPIC -L$(CUDD_DIR)/dddmp -L$(CUDD_DIR)/epd -L$(CUDD_DIR)/util -L$(CUDD_DIR)/st -L$(CUDD_DIR)/mtr -L$(CUDD_DIR)/cudd -L$(CUDD_DIR)/include
  INCLUDE    =  -I$(CUDD_DIR)/include  
  LIBS       = -lcudd -lmtr -lst -lutil -lepd -ldddmp

all: system.o pessoa2.out

system.o: System.cpp System.h
	 		$(CXX) $(CFLAGS) $(INCLUDE) $(LIBS) -stdlib=libc++ -c System.cpp

pessoa2.out: Pessoa2.cpp  
	 		$(CXX) $(CFLAGS) $(INCLUDE) system.o $(LIBS) -stdlib=libc++ Pessoa2.cpp -o ../pessoa2.out
	 		
